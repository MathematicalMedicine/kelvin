#!perl -w
use strict;
use KelvinConfig;
use KelvinDataset;
use KelvinFamily;
use DBI;

my $batchsize = 17;
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
my $pedfile;
my $family;
my $markercount;
my $serverid = undef;
my $models;
my $key;

my %positions;
my %reversemap;
my $merlinmodelid;
my $pedid;
my $position;
my @parts;
my @args;
my $merlinfailed = 0;

my %markerLK_modelids;
my @modelids;
my @LKs;

my $dsn;
my $dbh;
my $aref;
my $href;
my $line;

print ("$0 starting on $ENV{HOSTNAME} in $ENV{PWD}, pid $$". (exists ($ENV{JOB_ID}) ? ", job ID $ENV{JOB_ID}" : ""). "\n");

(exists ($ENV{JOB_ID})) and $prefix = $ENV{JOB_ID};
(exists ($ENV{SGE_TASK_ID}) && $ENV{SGE_TASK_ID} ne 'undefined')
    and $prefix .= '.' . $ENV{SGE_TASK_ID};

($configfile && -f $configfile)
    or die ("usage: $0 <configfile>\n");
print ("Using configfile $configfile\n");
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
print ("StudyId $study{id}");
print (($liability > 1) ? ", 1 liability class" : ", $liability liability classes");
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
# Release the KelvinDataset object
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
    $predataset = undef;
}
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
print ("ServerId is $serverid\n");

while (1) {
    $models = db_get_model_batch (\%study, $serverid, \%markerLK_modelids);
    (scalar (keys (%$models)) == 0) and last;
    print ("Selected ". scalar (keys (%$models)). " models\n");

    # $models{modelDetails}{pedigreeID}{"trait"|position} = ModelId
    # modelDetails := dgf-pen1dd-pen1Dd[-pen1dD]-pen1DD[-pen2dd...]

    %positions = %reversemap = @args = ();
    $merlinmodelid = 1;
    open (FH, ">$prefix.models") or die ("open '$prefix.models' failed, $!\n");
    foreach $key (keys (%$models)) {
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

    print ("Running '". join (' ', @args). "'\n");
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
	}
	if (exists ($$models{$key}{$parts[1]}{$parts[2]})) {
	    push (@modelids, $$models{$key}{$parts[1]}{$parts[2]});
	    push (@LKs, $parts[6]);
	    #print (OUT "combined $parts[1] $parts[2] $key $parts[6]\n");
	    delete ($$models{$key}{$parts[1]}{$parts[2]});
	}

	if (exists ($markerLK_modelids{$parts[1]}{$parts[2]})) {
	    push (@modelids, $markerLK_modelids{$parts[1]}{$parts[2]});
	    push (@LKs, $parts[5]);
	    #print (OUT "marker $parts[1] $parts[2] $parts[5]\n");
	    delete ($markerLK_modelids{$parts[1]}{$parts[2]});
	}
    }
    close (FH);
    #close (OUT);

    foreach $key (keys (%$models)) {
	foreach $pedid (keys(%{$$models{$key}})) {
	    foreach (keys (%{$$models{$key}{$pedid}})) {
		print ("leftover model: $$models{$key}{$pedid}{$_} $pedid $_ $key\n");
	    }
	}
    }

    foreach $pedid (keys (%markerLK_modelids)) {
	foreach $position (keys (%{$markerLK_modelids{$pedid}})) {
	    print ("leftover pedid in markerLK_modelids: $pedid at pos $position\n");
	}
    }

    db_update_LKs (\%study, \@modelids, \@LKs, 100);
#    unlink ("$prefix.models", "$prefix.par");

    my $va = 1;
    while (-f "$prefix.models.$va") {
	$va++;
    rename ("$prefix.models", "$prefix.models.$va");
    rename ("$prefix.par", "$prefix.par.$va");
    }
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
    
    (($sth = $dbh->prepare ("select distinct (PedigreeSId) ".
			    "from PedigreePositions ".
			    "where StudyId = ? ".
			    "  and PedigreeSId regexp ? ".
			    "  and PedigreeSId not regexp ? ".
			    "order by PedigreeSId limit 1")) &&
     $sth->execute (@$study{qw/id pedregex pednotregex/}))
	or die ("DBI select from PedigreePositions failed, $DBI::errstr\n");
    $aref = $sth->fetchrow_arrayref;
    $$study{minPedigreeSId} = $$aref[0];
    $sth->finish;

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
    my ($study, $serverid, $markerLK_modelids) = @_;
    my $sth;
    my $partcols = ", d1.DGF";
    my $parttabs = "";
    my $whereclause1 = "";
    my $whereclause2 = " and d1.DGF != -1";
    my @mpidcols = ();
    my @mpids = ();
    my @modelcols = (3);
    my @modelids = ();
    my @tmpmodelids = ();
    my @status;
    my $href = {};
    my $aref;
    my $rowref;
    my $offset;
    my $key;

    db_connect ($study);

    # Build some where clauses for selecting trait models given the number of liability
    # classes in the analysis, and whether or not imprinting is enabled. NOTE that we
    # select penetrances out of the database in the  order in which Merlin will want to
    # see them, which is mostly the reverse of how Kelvin would want them.

    if ($$study{ImprintingFlag} !~ /^y$/i) {
	map {
	    $offset = ($_ - 1) * 4;
	    $partcols .= ", d$_.MPId, d$_.LittlePen, d$_.BigLittlePen, d$_.BigPen";
	    push (@mpidcols, 4 + $offset);
	    push (@modelcols, 5 + $offset, 6 + $offset, 7 + $offset);
	} (1 .. $$study{LiabilityClassCnt});
    } else {
	map {
	    $offset = ($_ - 1) * 5;
	    $partcols .= ", d$_.MPId, d$_.LittlePen, d$_.BigLittlePen, d$_.LittleBigPen, d$_.BigPen";
	    push (@mpidcols, 4 + $offset);
	    push (@modelcols, 5 + $offset, 6 + $offset, 7 + $offset, 8 + $offset);
	} (1 .. $$study{LiabilityClassCnt});
    }
    map {
	$parttabs .= ", DModelParts d$_";
	$whereclause1 .= " and m.LC${_}MPId = d$_.MPId ";
	$whereclause2 .= " and m.LC${_}MPId = ? ";
    } (1 .. $$study{LiabilityClassCnt});

    # Any given trait model will apply to every pedigree (that is, there will be ModelIds
    # for that trait model for every pedigree). We only select ModelIds for one pedigree
    # (the lowest numbered pedigree that matches our configured regexs) to limit the
    # number of rows we need to lock in the database. We'll use trait model details to
    # select ModelIds for these models and all pedigrees later on.

    ($sth = $dbh->prepare ("select m.ModelId, p.PedigreeSId, ".
			    "  p.PedTraitPosCM". $partcols. " ".
			    "from Models m, PedigreePositions p". $parttabs. " ".
			    "where p.StudyId = ? ".
			    "  and p.PedigreeSId = ? ".
			    "  and p.PedTraitPosCM = -9999.99 ".
			    "  and m.ServerId is NULL ".
			    "  and m.Likelihood is NULL ".
			    "  and m.PedPosId = p.PedPosId ". $whereclause1.
			    "limit $batchsize for update"))
	or die ("DBI prepare select trait LK models for pedigree $study{minPedigreeSId} failed, $DBI::errstr\nSQL is '". $dbh->{Statement}. "'\n");
    #print ("SQL is ", $sth->{Statement}, "\n");
    while (1) {
	if (! $sth->execute (@$study{qw/id minPedigreeSId/})) {
	    die ("DBI select trait LK models for pedigree $study{minPedigreeSId} failed, $DBI::errstr\nSQL is '". $dbh->{Statement}. "'\n");
	} elsif (defined ($DBI::errstr) && $DBI::errstr =~ /try restarting transaction/) {
	    print ("retrying select trait LK models for pedigree $study{minPedigreeSId}\n");
	} else {
	    last;
	}
    }
    $aref = $sth->fetchall_arrayref;
    $sth->finish;

    # We get ModelIds for each trait model for the 'first family', and store those in
    # a list so we can update ownership, below. We also get MPIds (Model Part Ids) for
    # the models, so we can select the same trait models for all the other pedigrees
    # later. And we store each ModelId in a hash indexed by model details (DGF and
    # penetrance vector) and pedigree.

    my %pedids = ();

    foreach $rowref (@$aref) {
	push (@modelids, $$rowref[0]);
	push (@mpids, [ @$rowref[@mpidcols] ]);
	$key = join ('-', @$rowref[@modelcols]);
	$$href{$key}{$$rowref[1]}{trait} = $$rowref[0];
	$pedids{$$rowref[1]} = '';
    }

    print ("first family trait model count: ". scalar (@modelids). ", total keys: ". scalar (keys (%$href)). ", total pedids: ". scalar (keys (%pedids))."\n");

    # Update the current batch of models with our ServerId
    db_update_modelids ($study, $serverid, \@modelids);

    # Now select the markers-only ModelIds. There will be one of these for each pedigree
    # and calculation position. Since LG uses all markers at once, each pedigree will
    # have the same markers-only likelihood at every position. We store each ModelId in
    # a list for updating model ownership, below. We also store each ModelId in a hash
    # indexed by pedigree.

    %$markerLK_modelids = ();
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
    $aref = $sth->fetchall_arrayref;
    $sth->finish;

    %pedids = ();

    foreach $rowref (@$aref) {
	$$markerLK_modelids{$$rowref[1]}{$$rowref[2]} = $$rowref[0];
	push (@modelids, $$rowref[0]);
	exists ($pedids{$$rowref[1]}) or $pedids{$$rowref[1]} = 0;
	$pedids{$$rowref[1]}++;
    }

    if (scalar (@modelids)) {
	my $count = $pedids{(keys (%pedids))[0]};
	if ($count * scalar (keys (%pedids)) != scalar (@modelids)) {
	    print ("not all families returned the same number of marker models:\n");
	    print (map { "$_ => $pedids{$_}\n" } keys (%pedids));
	}
	print ("marker model count: ". scalar (@modelids). "\n");

	# Update the current batch of models with our ServerId
	db_update_modelids ($study, $serverid, \@modelids);
    } else {
	print ("marker model count: 0\n");
    }

    # Now select all the trait-only and combined ModelIds for all the families, based on
    # the Model Part Ids we selected above. Leave out the trait-only models for the 'first
    # family' that we've already selected.

    $sth = $dbh->prepare ("select m.ModelId, p.PedigreeSId, ".
			  "  p.PedTraitPosCM". $partcols. " ".
			  "from Models m, PedigreePositions p". $parttabs. " ".
			  "where p.StudyId = ? ".
			  "  and p.PedigreeSId regexp ? ".
			  "  and p.PedigreeSId not regexp ? ".
			  "  and not (p.PedTraitPosCM = -9999.99 ".
			  "           and p.PedigreeSId = ?) ".
			  "  and m.ServerId is NULL ".
			  "  and m.Likelihood is NULL ".
			  "  and m.PedPosId = p.PedPosId ".
			  $whereclause1 . $whereclause2)
	or die ("DBI prepare select trait/combined LK models failed, $DBI::errstr\nSQL is '". $dbh->{Statement}. "'\n");
    #print ("SQL is ", $sth->{Statement}, "\n");

    foreach (@mpids) {
	$sth->execute (@$study{qw/id pedregex pednotregex minPedigreeSId/}, @$_)
	    or die ("DBI execute select trait/combined LK models failed, $DBI::errstr\n");
	
	$aref = $sth->fetchall_arrayref;
	(defined ($aref) && scalar (@$aref) > 0)
	    or die ("select model parts returned no rows for MPIds ". join (',', $@). "\n");
	foreach $rowref (@$aref) {
	    push (@modelids, $$rowref[0]);
	    $key = join ('-', @$rowref[@modelcols]);
	    if ($$rowref[2] eq '-9999.99') {
		$$href{$key}{$$rowref[1]}{trait} = $$rowref[0];
	    } else {
		$$href{$key}{$$rowref[1]}{$$rowref[2]} = $$rowref[0];
	    }
	}
    }
    $sth->finish;

    print ("remaining model count: ", scalar (@modelids), "\n");
    
    # Update the current batch of models with our ServerId
    db_update_modelids ($study, $serverid, \@modelids, 100);

    return ($href);
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
    print ("updated $count model IDs with $retries retries\n");
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
