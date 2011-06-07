#!perl -w
use strict;
use POSIX qw(ceil);
use KelvinConfig;
use KelvinDataset;
use KelvinFamily;
use DBI;

my $taskcount = 9;
my $batchsize;
my $prefix = 'merlin';

my $configfile = shift (@ARGV);
my $config;
my %study;
my $imprinting = 0;
my $liability;
my $traitname = undef;
my $lcname = undef;

my $dataset;
my $predataset;
my $family;
my $markercount;
my $serverid = undef;
my $models;
my $key;

my %positions;
my %reversemap;
my $merlinmodelid;
my $pedid;
my @parts;
my @args;
my $merlinfailed = 0;

my @modelids;
my @LKs;
my ($trait_cnt, $combined_cnt, $marker_cnt) = (0, 0, 0);

my $dsn;
my $dbh;
my $aref;
my $href;
my $line;

print (ts(), "$0 starting on $ENV{HOSTNAME} in $ENV{PWD}, pid $$". (exists ($ENV{JOB_ID}) ? ", job ID $ENV{JOB_ID}" : ""). "\n");

(exists ($ENV{JOB_ID})) and $prefix = $ENV{JOB_ID};
(exists ($ENV{SGE_TASK_ID}) && $ENV{SGE_TASK_ID} ne 'undefined')
    and $prefix .= '.' . $ENV{SGE_TASK_ID};
(exists ($ENV{SGE_TASK_LAST}) && $ENV{SGE_TASK_LAST} ne 'undefined')
    and $taskcount = $ENV{SGE_TASK_LAST};

# $taskcount is initialized to 9, up top, and is overridden if SGE_TASK_LAST is set
$batchsize = ceil (153 / $taskcount);

($configfile && -f $configfile)
    or die ("usage: $0 <configfile>\n");
print (ts(), "Using configfile $configfile\n");
$config = KelvinConfig->new ($configfile)
    or die ("KelvinConfig->new failed: $KelvinConfig::errstr\n");
$aref = $config->isConfigured ("Study")
    or die ("no Study directive in config file '$configfile'\n'");
@study{qw/id role dbhost dbname dbuser dbpasswd pedregex pednotregex/} = split (' ', $$aref[0]);
($study{role} =~ /server/i)
    or die ("Study role in config file '$configfile' is not 'server'\n");

$aref = $config->isConfigured ("Imprinting");
defined ($aref) and $imprinting = 1;
$aref = $config->isConfigured ("LiabilityClasses");
$liability = defined ($aref) ? $$aref[0] : 1;
print (ts(), "StudyId $study{id}");
print (($liability > 1) ? ", $liability liability classes" : ", 1 liability class");
print (($imprinting) ? ", imprinting enabled\n" : "\n");

$$href{LocusFile} = $ {$config->isConfigured ("LocusFile")}[0];
$$href{PedigreeFile} = $ {$config->isConfigured ("PedigreeFile")}[0];
$dataset = KelvinDataset->new ($href)
    or die ("new KelvinDataset failed, $KelvinDataset::errstr\n");
$study{markercount} = scalar (@{$dataset->markerOrder});
($study{markercount} > 0) or die ("no markers in dataset\n");
foreach (@{$dataset->traitOrder}) {
    $href = $dataset->getTrait ($_);
    ($$href{flag} eq 'A' && ! defined ($traitname)) and $traitname = $_;
    ($$href{flag} eq 'C' && ! defined ($lcname)) and $lcname = $_;
}
(defined ($traitname)) or die ("no affection status in dataset\n");
(defined ($lcname) || $liability == 1) or die ("no liability class in dataset\n");

# If we've been given a post-MAKEPED pedigree file, we need to create a pre-makeped file
# for Merlin. Also, we'll create a new pedigree file if the regular expressions in the Study
# directive actually exclude any pedigrees, so we only give Merlin the pedigrees we want
# to analyze. NOTE this is actually a bad design. Preprocessing ped files should be done
# before this script is called and it's not clear that InitStudy properly honors pedigree
# exclusion via regex anyway.

($family = $dataset->readFamily)
    or die ("KelvinDataset readFamily failed: $KelvinDataset::errstr\n");
if ($study{pedregex} ne "^.*\$" || $study{pednotregex} ne 'xyzzy' || $family->origfmt ne 'pre') {
    $config->setDirective ("PedigreeFile", $ {$config->isConfigured ("PedigreeFile")}[0].'.pre');
    $predataset = $dataset->copy;
    $predataset->writePedigreefile ({pedigreefile => $ {$config->isConfigured ("PedigreeFile")}[0],
				     premakeped => 1})
	or die ("KelvinDataset writePedigreefile failed, $KelvinDataset::errstr\n");
    while ($family) {
	if ($family->pedid =~ /$study{pedregex}/ && $family->pedid !~ /$study{pednotregex}/) {
	    $family = $family->map ($predataset);
	    $family->write;
	}
	$family = $dataset->readFamily;
    }
    (defined ($family)) 
	or die ("KelvinDataset readFamily failed: $KelvinDataset::errstr\n");
    # Release the pre-makeped copy of the KelvinDataset object to free memory
    $predataset = undef;
}
# Release the KelvinDataset object to free memory
$dataset = undef;

$dsn = "DBI:mysql:host=$study{dbhost};database=$study{dbname}";

db_get_study_info (\%study);
(($study{ImprintingFlag} =~ /^y$/i) == ($imprinting == 1))
    or die ("Imprinting flag mismatch between config and database\n");
($study{LiabilityClassCnt} == $liability)
    or die ("Liability Class mismatch between config and database\n");
(defined ($study{PendingWorkFlag}) && $study{PendingWorkFlag} eq 'D')
    and die ("Study PendingWorkFlag is already set to 'done'\n");
$serverid = db_get_serverid (\%study);
print (ts(), "ServerId is $serverid\n");

while (1) {
    $models = db_get_model_batch (\%study, $serverid);
    (scalar (keys (%$models)) == 0) and last;
    print (ts(), "Selected ". scalar (keys (%$models)). " models\n");

    # $models{modelDetails}{pedigreeID}{"trait"|position} = ModelId  -- or --
    # $models{"marker"}{pedigreeID}{position} = ModelId 
    # modelDetails := dgf-pen1dd-pen1Dd[-pen1dD]-pen1DD[-pen2dd...]

    %positions = %reversemap = @args = ();
    $merlinmodelid = 1;
    open (FH, ">$prefix.models") or die ("open '$prefix.models' failed, $!\n");
    foreach $key (keys (%$models)) {
	($key eq 'marker') and next;
	@parts = split (/-/, $key);
	if ($liability == 1) {
	    print (FH join (' ', $traitname, shift (@parts), merlin_pens (\@parts, $imprinting), "model-$merlinmodelid"), "\n");
	} else {
	    print (FH join (' ', $traitname, shift (@parts), '*', "model-$merlinmodelid"), "\n");
	    map {
		print (FH "$lcname = $_ ". merlin_pens (\@parts, $imprinting). "\n");
	    } (1 .. $liability - 1);
	    print (FH "OTHERWISE ". merlin_pens (\@parts, $imprinting). "\n");
	}
	$reversemap{"model-".$merlinmodelid++} = $key;

	# We don't really care about the values here, we only need to make sure that
	# the keys of %positions is a superset of all the positions that were selected
	# for any model/family combination.
	foreach $pedid (keys (%{$$models{$key}})) {
	    @positions{keys (%{$$models{$key}{$pedid}})} = '';
	}
    }
    close (FH);
    (exists ($positions{trait})) and delete ($positions{trait});

    push (@args, ($config->isConfigured ("SexLinked") ? "minx" : "merlin"));
    push (@args, '-d', $ {$config->isConfigured ("LocusFile")}[0]);
    push (@args, '-f', $ {$config->isConfigured ("FrequencyFile")}[0]);
    push (@args, '-m', $ {$config->isConfigured ("MapFile")}[0]);
    push (@args, '-p', $ {$config->isConfigured ("PedigreeFile")}[0]);
    push (@args, '--model', "$prefix.models");
    push (@args, '--prefix', $prefix);
    ($imprinting) and push (@args, "--parentOfOrigin");
    push (@args, '--bits', '80', '--quiet', '--perFamily', '--likelihoodServer');
    ($imprinting) and push (@args, '--parentOfOrigin');
    push (@args, '--positions', join (',', sort ({$a <=> $b} keys (%positions))));

    print (ts(), "Running '". join (' ', @args). "'\n");
    open (FH, join (' ', @args, '2>&1 |'))
	or die ("open '$args[0]' failed, $!\n");
    while ($line = <FH>) {
	($line =~ /(fatal|crash|skipped)/i) and $merlinfailed = 1;
	print ($line);
    }
    close (FH)
	or die ("close '$args[0]' failed, $!\n");
    ($? != 0) and die ("args[0] exited with error status ", $? >> 8, "\n");
    ($merlinfailed) and die ("$args[0] reported a non-fatal failure\n");

    @modelids = @LKs = ();
    open (FH, "$prefix.par") or die ("open '$prefix.par' failed, $!\n");
    #open (OUT, ">$prefix.debug") or die ("open '$prefix.debug' failed, $!\n");
    while ($line = <FH>) {
	$line =~ /^MODEL/ and next;
	@parts = split (' ', $line);
	# fields are: modelname pedigreeId position LOD traitLK markerLK combinedLK
	
	$key = $reversemap{$parts[0]};
	if (exists ($$models{$key}{$parts[1]}{trait})) {
	    push (@modelids, $$models{$key}{$parts[1]}{trait});
	    push (@LKs, $parts[4]);
	    #print (OUT "trait $parts[1] $key $parts[4]\n");
	    delete ($$models{$key}{$parts[1]}{trait});
	    $trait_cnt++;
	}
	if (exists ($$models{$key}{$parts[1]}{$parts[2]})) {
	    push (@modelids, $$models{$key}{$parts[1]}{$parts[2]});
	    push (@LKs, $parts[6]);
	    #print (OUT "combined $parts[1] $parts[2] $key $parts[6]\n");
	    delete ($$models{$key}{$parts[1]}{$parts[2]});
	    $combined_cnt++;
	}
	if (exists ($$models{marker}{$parts[1]}{$parts[2]})) {
	    push (@modelids, $$models{marker}{$parts[1]}{$parts[2]});
	    push (@LKs, $parts[5]);
	    #print (OUT "marker $parts[1] $parts[2] $key $parts[5]\n");
	    delete ($$models{marker}{$parts[1]}{$parts[2]});
	    $marker_cnt++;
	}
    }
    close (FH);
    #close (OUT);
    print (ts(), "found LKs for $trait_cnt trait models, $combined_cnt cobined models, and $marker_cnt marker models\n");

    foreach $key (keys (%$models)) {
	foreach $pedid (keys(%{$$models{$key}})) {
	    foreach (keys (%{$$models{$key}{$pedid}})) {
		print (ts(), "leftover model: $$models{$key}{$pedid}{$_} $pedid $_ $key\n");
	    }
	}
    }

    db_update_LKs (\%study, \@modelids, \@LKs, 100);
    unlink ("$prefix.models", "$prefix.par");

#    my $va = 1;
#    while (-f "$prefix.models.$va") {
#	$va++;
#	rename ("$prefix.models", "$prefix.models.$va");
#	rename ("$prefix.par", "$prefix.par.$va");
#    }
}

END {
    my $status = $?;

    if (defined ($serverid)) {
	db_set_server_status (\%study, $serverid, $status);
    }

    $dbh->disconnect;
    exit ($status);
}


sub merlin_pens
{
    my ($parts, $imprinting) = @_;
    
    return (join (',', splice (@$parts, 0, ($imprinting ? 4 : 3))));
}


sub db_connect
{
    my ($study) = @_;

    $dbh = DBI->connect_cached ($dsn, $$study{dbuser}, $$study{dbpasswd},
#				{AutoCommit => 0})
				{AutoCommit => 0, PrintError => 0})
	or die ("DBI connect to '$dsn' as $study{dbuser} failed, $DBI::errstr\n");
}


sub db_get_study_info
{
    my ($study) = @_;
    my $sth;
    my $href;
    my $aref;
    
    db_connect ($study);
    (($sth = $dbh->prepare ("select LiabilityClassCnt, ImprintingFlag, PendingWorkFlag ".
			    "from Studies where StudyId = ?")) &&
     $sth->execute ($$study{id}))
	or die ("DBI select from Studies failed, $DBI::errstr\n");
    $href = $sth->fetchrow_hashref;
    $sth->finish;
    @$study{keys (%$href)} = @$href{keys (%$href)};
    
    return (1);
}


sub db_get_serverid
{
    my ($study) = @_;
    my $sth;
    my $aref;

    db_connect ($study);
    (($sth = $dbh->prepare ("insert into Servers ".
			    "  (ServerId, Hostname, ProcessId, StudyId, ".
			    "   PedigreeRegEx, PedigreeNotRegEx, MarkerCount, ".
			    "   Algorithm, ProgramVersion, ListenerSocketId) ".
			    "values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")) &&
     $sth->execute (0, $ENV{HOSTNAME}, $$, @$study{qw/id pedregex pednotregex markercount/},
		    "LG", "1.1.2", undef))
	or die ("DBI insert into Servers failed, $DBI::errstr\n");
    $dbh->commit
	or die ("DBI commit insert into Servers failed, $DBI::errstr\n");
    $aref = $dbh->selectall_arrayref ("select LAST_INSERT_ID()")
	or die ("DBI select ServerID failed, $DBI::errstr\n");
    (defined ($$aref[0]) && defined ($$aref[0][0]) && $$aref[0][0] =~ /^\d+$/)
	or die ("ServerId returned from DB is garbage\n");
    $serverid = $$aref[0][0];
    $sth->finish;
    return ($serverid);
}



sub db_get_model_batch
{
    my ($study, $serverid) = @_;
    my $batch = {};

    my $MPId_colnames = '';
    my $DModelPart_colnames = '';
    my $DModelPart_tabnames = '';
    my $MPId_whereclause = '';

    my @MPIds = ();
    my @Models = ();
    my @modelids = ();
    my ($modelid, $pedid, $traitpos, $key);
    my $sth;
    my $aref;

    # This is a cheap hack: with really fast datasets, sometimes the higher-numbered
    # (n > 1) tasks snorkle up all the trait-marker models, and the first task gets
    # none. This is a problem, since only the first task will select the marker-only
    # models, but only if gets some trait-marker models also. So, we give the first
    # task a five second head start.
    (defined ($ENV{SGE_TASK_ID}) && $ENV{SGE_TASK_ID} != 1)
	and sleep (5);

    $MPId_colnames = join (', ', map { "LC${_}MPId" } (1 .. $$study{LiabilityClassCnt}));
    
    print (ts(), "Selecting model part IDs\n");
    (($sth = $dbh->prepare ("select LGModelId, $MPId_colnames from LGModels ".
			    "where StudyId = ? and ServerId is NULL ".
			    "limit $batchsize for update"))
     && $sth->execute ($study{id}))
	or die ("select from LGModels failed, $DBI::errstr\n");
    # @MPIds will contain a list of arrayrefs. Each arrayref contains one or more model
    # part Ids, depending on the number of liability classes
    @MPIds = @{$sth->fetchall_arrayref};
    $sth->finish;
    print (ts(), "selected ". scalar (@MPIds). " rows\n");
    (scalar (@MPIds) == 0) and return ($batch);

    foreach (@MPIds) {
	push (@modelids, shift (@$_));
    }
    db_update_lgmodels ($study, $serverid, \@modelids);
    @modelids = ();

    # select trait/combined modelIDs from Models based on the LCxMPIds we got from LGModels

    map {
	$DModelPart_tabnames .= ", DModelParts d$_";
	$MPId_whereclause .= " and m.LC${_}MPId = d$_.MPId and m.LC${_}MPId = ?";
    } (1 .. $$study{LiabilityClassCnt});

    if ($$study{ImprintingFlag} !~ /^y$/i) {
	map {
	    $DModelPart_colnames .= ", d$_.LittlePen, d$_.BigLittlePen, d$_.BigPen";
	} (1 .. $$study{LiabilityClassCnt});
    } else {
	map {
	    $DModelPart_colnames .= ", d$_.LittlePen, d$_.BigLittlePen, d$_.LittleBigPen, d$_.BigPen";
	} (1 .. $$study{LiabilityClassCnt});
    }

    # We're guaranteed to have at least 'DModelParts d1' in $DModelPart_tabnames, so it's
    # safe to refer to 'd1.DGF' in the static part of the SQL. Also, merlin_lk_prepare.pl
    # already excluded marker-only LK models (that is, DGF == -1), so we don't need to 
    # exclude those again here.

    print (ts(), "Selecting trait LK and combined LK models\n");
    ($sth = $dbh->prepare ("select m.ModelId, p.PedigreeSId, p.PedTraitPosCM, d1.DGF".
			   $DModelPart_colnames . " ".
			   "from Models m, PedigreePositions p".
			   $DModelPart_tabnames . " ".
			   "where p.StudyId = ? ".
			   "  and p.PedigreeSId regexp ? ".
			   "  and p.PedigreeSId not regexp ? ".
			   "  and p.PedPosId = m.PedPosId".
			   $MPId_whereclause))
	or die ("prepare select from Models failed, $DBI::errstr\n");
    print ("SQL is ". $sth->{Statement}. "\n");
    foreach (@MPIds) {
	print ("Selecting for model part IDs ". join (', ', @$_). "\n");
	# Recall that each element in @MPIds is an arrayref containing one or more model part IDs
	$sth->execute (@$study{qw/id pedregex pednotregex/}, @$_)
	    or die ("execute select from Models failed, $DBI::errstr\n");
	@Models = @{$sth->fetchall_arrayref};
	(scalar (@Models) == 0) and die ("select returned no models\n");
	foreach $aref (@Models) {
	    ($modelid, $pedid, $traitpos) = splice (@$aref, 0, 3);
	    $key = join ('-', @$aref);
	    if ($traitpos == -9999.99) {
		$$batch{$key}{$pedid}{trait} = $modelid;
	    } else {
		$$batch{$key}{$pedid}{$traitpos} = $modelid;
	    }
	    push (@modelids, $modelid);
	}
	print ("Select returned ". scalar (@Models). " rows\n");
    }
    $sth->finish;
    db_update_modelids ($study, $serverid, \@modelids);
    @modelids = ();

    # This next select is a great opportunity for database deadlocking, so we're only
    # going to try it if 1) we're apparently running as a standalone (not an array) job,
    # or 2) if we are part of an array job, but we're running as task ID 1.

    if (defined ($ENV{SGE_TASK_ID}) && $ENV{SGE_TASK_ID} != 1) {
	print (ts(), "Skipping marker LK models\n");
	return ($batch);
    }

    print (ts(), "Selecting marker LK models\n");
    ($sth = $dbh->prepare ("select m.ModelId, p.PedigreeSId, p.PedTraitPosCM ".
			   "from Models m, PedigreePositions p, DModelParts d ".
			   "where p.StudyId = ? ".
			   "  and p.PedigreeSId regexp ? ".
			   "  and p.PedigreeSId not regexp ? ".
			   "  and p.PedPosId = m.PedPosId ".
			   "  and m.ServerId is NULL ".
			   "  and m.Likelihood is NULL ".
			   "  and m.LC1MPId = d.MPId ".
			   "  and d.DGF = -1 ".
			   "for update"))
	or die ("DBI prepare select marker LK models failed, $DBI::errstr\n");
    while (1) {
	if (! $sth->execute (@$study{qw/id pedregex pednotregex/})) {
	    die ("DBI select marker LK models failed, $DBI::errstr\n");
	} elsif (defined ($DBI::errstr) && $DBI::errstr =~ /try restarting transaction/) {
	    print ("retrying select marker LK models\n");
	} else {
	    last;
	}
    }
    @Models = @{$sth->fetchall_arrayref};
    $sth->finish;
    print (ts(), "select returned ". scalar (@Models). " rows\n");
    (scalar (@Models) == 0) and return ($batch);

    foreach $aref (@Models) {
        ($modelid, $pedid, $traitpos) = @$aref;
	$$batch{marker}{$pedid}{$traitpos} = $modelid;
	push (@modelids, $modelid);
    }
    db_update_modelids ($study, $serverid, \@modelids);
    @modelids = ();
    return ($batch);
}


sub db_update_lgmodels
{
    my ($study, $serverid, $LGModelIds) = @_;
    my $sth;
    my @status;
    my @errors;
    my $retries = 0;

    ($sth = $dbh->prepare ("update LGModels set ServerId = $serverid ".
			   "where LGModelId = ?"))
	or die ("prepare update LGModels failed, $DBI::errstr\n");
    while (1) {
	if ($sth->execute_array ({ArrayTupStatus => \@status}, $LGModelIds)) {
	    last;
	} elsif (! all_retry_errors (\@status, \@errors)) {
	    die ("execute update Models failed:\n", map { "$_\n" } @errors);
	}
	$retries++;
    }
    $dbh->commit or die ("commit update LGModels failed, $DBI::errstr\n");
    $sth->finish;
    print (ts(), "Claimed ". scalar (@$LGModelIds). " LGModels with $retries retries\n");
    return (1);
}


sub db_update_modelids
{
    my ($study, $serverid, $modelids, $batchsize) = @_;
    my $sth;
    my @status;
    my @batch;
    my @errors;
    my $count = 0;
    my $retries = 0;

    db_connect ($study);
    ($sth = $dbh->prepare ("update Models set ServerId = $serverid, ".
			   "  MarkerCount = $$study{markercount}, ".
			   "  StartTime = NOW() ".
			   "where ModelId = ?"))
	or die ("DBI prepare update Models failed, $DBI::errstr\n");

    (defined ($batchsize)) or $batchsize = scalar (@$modelids);
    while (scalar (@$modelids)) {
	@batch = splice (@$modelids, 0, $batchsize);
	while (1) {
	    if ($sth->execute_array ({ArrayTupleStatus => \@status}, \@batch)) {
		#print ("updated ", scalar (@batch), " Models\n");
		$count += scalar (@batch);
		last;
	    } elsif (! all_retry_errors (\@status, \@errors)) {
		die ("execute update Models failed:\n", map { "$_\n" } @errors);
	    }
	    #print ("retrying transaction\n");
	    $retries++;
	}
	$dbh->commit or die ("commit update Models failed, $DBI::errstr\n");
    }
    $sth->finish;
    print (ts(), "Updated $count model IDs with $retries retries\n");
    return (1);
}


sub db_update_LKs
{
    my ($study, $modelids, $LKs, $batchsize) = @_;
    my $sth;
    my @status;
    my @idbatch;
    my @LKbatch;
    my @errors;
    my $count = 0;
    my $retries = 0;

    db_connect ($study);
    ($sth = $dbh->prepare ("update Models ".
			   "set Likelihood = ?, EndTime = NOW() ".
			   "where ModelId = ?"))
	or die ("DBI prepare update of model LKs failed, $DBI::errstr\n");

    (defined ($batchsize)) or $batchsize = scalar (@$modelids);
    while (scalar (@$modelids)) {
	@idbatch = splice (@$modelids, 0, $batchsize);
	@LKbatch = splice (@$LKs, 0, $batchsize);
	while (1) {
	    if ($sth->execute_array ({ArrayTupleStatus => \@status}, \@LKbatch, \@idbatch)) {
		#print ("updated ", scalar (@idbatch), " model LKs\n");
		$count += scalar (@idbatch);
		last;
	    } elsif (! all_retry_errors (\@status, \@errors)) {
		die ("execute update model LKs failed:\n", map { "$_\n" } @errors);
	    }
	    #print ("retrying transaction\n");
	    $retries++;
	}
	$dbh->commit or die ("commit update model LKs failed, $DBI::errstr\n");
    }
    $sth->finish;
    print ("updated $count model LKs with $retries retries\n");
    return (1);
}


sub db_set_server_status
{
    my ($study, $serverid, $status) = @_;
    my $sth;

    db_connect ($study);
    (($sth = $dbh->prepare ("update Servers set StopTime = NOW(), ExitStatus = ? ".
			    "where ServerID = ?"))
     && $sth->execute ($status, $serverid))
	or warn ("DBI update Servers failed, $DBI::errstr\n");
    $dbh->commit or warn ("DBI commit update Servers failed, $DBI::errstr\n");
    return (1);
}


sub all_retry_errors
{
    my ($status, $errors) = @_;
    my $va;

    @$errors = ();
    for ($va = 0; $va < scalar (@$status); $va++) {
	(! defined ($$status[$va])) and $$status[$va] = [0, "skipped for unknown reasons"];
	ref ($$status[$va]) or next;
	$$status[$va][1] =~ /try restarting transaction/ and next;
	push (@$errors, $$status[$va][1]);
    }
    return (scalar (@$errors) ? undef : 1);
}

sub ts
{
    my @arr = localtime ();
    return (sprintf ("%02d/%02d %02d:%02d:%02d ", $arr[4]+1, @arr[3,2,1,0]));
}
