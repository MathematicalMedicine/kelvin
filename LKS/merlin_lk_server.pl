#!/usr/bin/env perl 
use strict;
use FindBin;
use lib split(/:+/, "!:$ENV{'TOOLPATH'}:$FindBin::Bin");
use warnings;
use POSIX qw(ceil);
use KelvinConfig;
use KelvinDataset;
use KelvinFamily;
use DBI;
use Data::Dumper;

my $taskcount = 9;
my $batchsize;
my $prefix = 'merlin';
my $svn_version='$Id$';

my $configfile = shift (@ARGV);
my $config;
my %study;
my $imprinting = 0;
my $liability;
my $traitname = undef;
my $lcname = undef;
my $traitflag;

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
my @modelLKs;
my @markersetids;
my @markersetLKs;
my ($trait_cnt, $combined_cnt, $marker_cnt) = (0, 0, 0);

my $dsn = undef;
my $dbh;
my $aref;
my $href;
my $line;
my $debug = 1;   # Used globally
my $auto = 1;    # Used globally

print (ts(), "$0 starting on $ENV{HOSTNAME} in $ENV{PWD}, pid $$". (exists ($ENV{JOB_ID}) ? ", job ID $ENV{JOB_ID}" : ""). "\n");
print (ts(), "Version $svn_version\n");
(exists($ENV{TOOLPATH})) and $ENV{PATH} = "$ENV{TOOLPATH}:$ENV{PATH}";
        # because PATH gets reset on some hosts for some reason, and we need it
        # properly set so we can invoke Merlin

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
@study{qw/label role dbhost dbname dbuser dbpasswd pedregex pednotregex/} = split (' ', $$aref[0]);
($study{role} =~ /server/i)
    or die ("Study role in config file '$configfile' is not 'server'\n");

db_get_study_info (\%study);

$aref = $config->isConfigured ("Imprinting");
defined ($aref) and $imprinting = 1;
$aref = $config->isConfigured ("LiabilityClasses");
$liability = defined ($aref) ? $$aref[0] : 1;
print (ts(), "StudyLabel '$study{label}', ID $study{id}");
print (($liability > 1) ? ", $liability liability classes" : ", 1 liability class");
print (($imprinting) ? ", imprinting enabled\n" : "\n");
if ($aref = $config->isConfigured ("QT")) {
    $study{TraitType} = 'QT';
    $traitflag = 'T';
} elsif ($aref = $config->isConfigured ("QTT")) {
    $study{TraitType} = 'QTT';
    $traitflag = 'T';
} else {
    $study{TraitType} = 'DT';
    $traitflag = 'A';
}
if (defined ($aref)) {
    if ($$aref[0] =~ /normal/i) {
	$study{TraitDist} = 'T';
    } elsif ($$aref[0] =~ /chisq/i) {
	$study{TraitDist} = 'ChiSq';
    } else {
	die ("unknown QT/QTT distribution '$$aref[0]' in config\n");
    }
}
print (ts(), "Trait is $study{TraitType}", exists ($study{TraitDist}) ? ", $study{TraitDist} distribution\n" : "\n");

(($study{ImprintingFlag} =~ /^y$/i) == ($imprinting == 1))
    or die ("Imprinting flag mismatch between config and database\n");
($study{LiabilityClassCnt} == $liability)
    or die ("Liability Class mismatch between config and database\n");
(defined ($study{PendingWorkFlag}) && $study{PendingWorkFlag} eq 'D')
    and die ("Study PendingWorkFlag is already set to 'done'\n");

$$href{LocusFile} = $ {$config->isConfigured ("LocusFile")}[0];
$$href{PedigreeFile} = $ {$config->isConfigured ("PedigreeFile")}[0];
$$href{MapFile} = $ {$config->isConfigured ("MapFile")}[0];
$dataset = KelvinDataset->new ($href)
    or die ("new KelvinDataset failed, $KelvinDataset::errstr\n");
$study{markercount} = scalar (@{$dataset->markerOrder});
$study{chromosome} = $dataset->chromosome;
($study{markercount} > 0) or die ("no markers in dataset\n");
foreach (@{$dataset->traitOrder}) {
    $href = $dataset->getTrait ($_);
    ($$href{flag} eq $traitflag && ! defined ($traitname)) and $traitname = $_;
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
#if ($study{pedregex} ne "^.*\$" || $study{pednotregex} ne 'xyzzy' || $family->origfmt ne 'pre') {
# Temporarily not considering pedregex and pednotregex
if ($family->origfmt ne 'pre') {
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

$serverid = db_get_serverid (\%study);
print (ts(), "ServerId is $serverid\n");

while (1) {
    $models = db_get_model_batch (\%study, $serverid);
    (scalar (keys (%$models)) == 0) and last;
    print (ts(), "Selected ". scalar (keys (%$models)). " models\n");

    # $models{modelDetails}{pedigreeID}{"trait"|position} = ModelId  -- or --
    # $models{"marker"}{pedigreeID}{position} = ModelId 
    # modelDetails (DT) = dgf:pen1dd:pen1Dd[:pen1dD]:pen1DD
    # modelDetails (QT-T) = dgf:SD1dd:SD1Dd[:SD1dD]:SD1DD:mean1dd:mean1Dd[:mean1dD]:mean1DD
    # modelDetails (QTT-T) = dgf:SD1dd:SD1Dd[:SD1dD]:SD1DD:mean1dd:mean1Dd[:mean1dD]:mean1DD:Threshold1
    # modelDetails (QT-ChiSq) = dgf:DoF1dd:DoF1Dd[:DoF1dD]:DoF1DD
    # modelDetails (QTT-ChiSq) = dgf:DoF1dd:DoF1Dd[:DoF1dD]:DoF1DD:Threshold1
    # for multiple liability classes, repeat everything except dgf for each class

    %positions = %reversemap = @args = ();
    $merlinmodelid = 1;
    open (FH, ">$prefix.models") or die ("open '$prefix.models' failed, $!\n");
    foreach $key (keys (%$models)) {
	($key eq 'marker') and next;
	@parts = split (/:/, $key);
	if ($liability == 1) {
	    print (FH join (' ', $traitname, shift (@parts), merlin_pens (\@parts, \%study), "model-$merlinmodelid"), "\n");
	} else {
	    print (FH join (' ', $traitname, shift (@parts), '*', "model-$merlinmodelid"), "\n");
	    map {
		print (FH "$lcname = $_ ". merlin_pens (\@parts, \%study). "\n");
	    } (1 .. $liability - 1);
	    print (FH "OTHERWISE ". merlin_pens (\@parts, \%study). "\n");
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
    ($study{TraitType} =~ /^QT/) and push (@args, "--qtDistribution $study{TraitDist}");
    ($study{TraitType} eq 'QTT') and push (@args, "--qtThreshold");
    push (@args, '--bits', '80', '--quiet', '--perFamily', '--likelihoodServer');
    push (@args, '--positions', join (',', sort ({$a <=> $b} keys (%positions))));
    (-f "merlin-clusters.freq") and push (@args, '--clusters', 'merlin-clusters.freq');

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

    @modelids = @modelLKs = @markersetids = @markersetLKs = ();
    open (FH, "$prefix.par") or die ("open '$prefix.par' failed, $!\n");
    #open (OUT, ">$prefix.debug") or die ("open '$prefix.debug' failed, $!\n");
    while ($line = <FH>) {
	$line =~ /^MODEL/ and next;
	@parts = split (' ', $line);
	# fields are: modelname pedigreeId position LOD traitLK markerLK combinedLK
	
	$key = $reversemap{$parts[0]};
	if (exists ($$models{$key}{$parts[1]}{trait})) {
	    push (@modelids, $$models{$key}{$parts[1]}{trait});
	    push (@modelLKs, $parts[4]);
	    #print (OUT "trait $modelids[-1] $parts[1] $key $parts[4]\n");
	    delete ($$models{$key}{$parts[1]}{trait});
	    $trait_cnt++;
	}
	if (exists ($$models{$key}{$parts[1]}{$parts[2]})) {
	    push (@modelids, $$models{$key}{$parts[1]}{$parts[2]});
	    push (@modelLKs, $parts[6]);
	    #print (OUT "combined $modelids[-1] $parts[1] $parts[2] $key $parts[6]\n");
	    delete ($$models{$key}{$parts[1]}{$parts[2]});
	    $combined_cnt++;
	}
	if (exists ($$models{marker}{$parts[1]}{$parts[2]})) {
	    push (@markersetids, $$models{marker}{$parts[1]}{$parts[2]});
	    push (@markersetLKs, $parts[5]);
	    #print (OUT "marker $markersetids[-1] $parts[1] $parts[2] $key $parts[5]\n");
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

    db_update_modelLKs (\%study, \@modelids, \@modelLKs, 100);
    db_update_markersetLKs (\%study, \@markersetids, \@markersetLKs, 100);
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
    my ($parts, $study) = @_;
    my $format = "%s,%s,%s";
    my $count = 3;

    if ($$study{ImprintingFlag} =~ /y/i) {
	$format .= ",%s";
	$count++;
    }

    if (exists ($$study{TraitDist}) && $$study{TraitDist} eq 'T') {
	$format = $format . " ". $format;
	$count += $count;
    }
    if ($$study{TraitType} eq 'QTT') {
	$format .= " %s";
	$count++;
    }
    return (sprintf ($format, splice (@$parts, 0, $count)));
}


sub db_connect
{
    my ($study) = @_;

    (defined ($dsn)) or $dsn = "DBI:mysql:host=$$study{dbhost};database=$$study{dbname}";

    $dbh = DBI->connect_cached ($dsn, $$study{dbuser}, $$study{dbpasswd},
				{AutoCommit => $auto, PrintError => 0})
	or die ("DBI connect to '$dsn' as $study{dbuser} failed, $DBI::errstr\n");
}


sub db_get_study_info
{
    my ($study) = @_;
    my $sth;
    my $href;
    my $aref;
    
    db_connect ($study);
    (($sth = $dbh->prepare ("select StudyId id, LiabilityClassCnt, ImprintingFlag, ".
			    "  PendingWorkFlag ".
			    "from Studies where StudyLabel = ?")) &&
     $sth->execute ($$study{label}))
	or die ("DBI select from Studies failed, $DBI::errstr\n");
    (defined ($href = $sth->fetchrow_hashref))
	or die ("DBI select from Studies returned no data\n");
    $sth->finish;
    @$study{keys (%$href)} = @$href{keys (%$href)};
    
    return (1);
}


sub db_get_serverid
{
    my ($study) = @_;
    my $serverid;
    my $sth;
    my $aref;

    db_connect ($study);
    (($sth = $dbh->prepare ("call ServerSignOn (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, \@sid)")) &&
     $sth->execute ($ENV{HOSTNAME}, $$, 1, @$study{qw/id pedregex pednotregex chromosome/},
		    "LG", $$study{markercount}, "1.1.2", 0, 0))
	or die ("DBI call ServerSignOn() failed, $DBI::errstr\n");
    $aref = $dbh->selectall_arrayref ("select \@sid")
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
    my $ModelPart_colnames = '';
    my $ModelPart_tabnames = '';
    my $MPId_whereclause = '';

    my @MPIds = ();
    my @Models = ();
    my @lgmodelids = ();
    my @modelids = ();
    my @markersetids = ();
    my ($modelid, $markersetid, $pedid, $traitpos, $key);
    my $colnames = '';
    my $sth;
    my $aref;

    # This is a cheap hack: with really fast datasets, sometimes the higher-numbered
    # (n > 1) tasks snorkle up all the trait-marker models, and the first task gets
    # none. This is a problem, since only the first task will select the marker-only
    # models, but only if gets some trait-marker models also. So, we give the first
    # task a five second head start.
    (defined ($ENV{SGE_TASK_ID}) && $ENV{SGE_TASK_ID} != 1)
	and sleep (5);

    $MPId_colnames = join (', ', map { "l.LC${_}MPId" } (1 .. $$study{LiabilityClassCnt}));
    
    print (ts(), "Selecting model part IDs\n");
    (($sth = $dbh->prepare ("select l.LGModelId, $MPId_colnames ".
                            "  from Analyses a, LGModels l ".
			    "where a.StudyId = ? ".
                            "  and a.PedigreeRegEx = ? ".
                            "  and a.PedigreeNotRegEx = ? ".
                            "  and l.StudyId = a.StudyId ".
                            "  and l.AnalysisId = a.AnalysisId ".
                            "  and ServerId is NULL ".
			    "limit $batchsize for update"))
     && $sth->execute (@$study{qw/id pedregex pednotregex/}))
	or die ("select from LGModels failed, $DBI::errstr\n");
    # @MPIds will contain a list of arrayrefs. Each arrayref contains one or more model
    # part Ids, depending on the number of liability classes
    @MPIds = @{$sth->fetchall_arrayref};
    $sth->finish;
    print (ts(), "Select returned ". scalar (@MPIds). " rows\n");
    (scalar (@MPIds) == 0) and return ($batch);

    foreach (@MPIds) {
        # This strips the LGModelId off of every row returned, leaving only
        # the model part IDs in each element of @MPIds.
	push (@lgmodelids, shift (@$_));
    }
    db_update_lgmodels ($study, $serverid, \@lgmodelids);

    # select trait/combined modelIDs from Models based on the LCxMPIds we got from LGModels

    if ($$study{TraitType} eq 'DT') {
	map {
	    $ModelPart_tabnames .= ", DModelParts d$_";
	} (1 .. $$study{LiabilityClassCnt});
	
	$colnames .= ", dXX.LittlePen, dXX.BigLittlePen";
	($$study{ImprintingFlag} =~ /^y$/i) and $colnames .= ", dXX.LittleBigPen";
	$colnames .= ", dXX.BigPen";

    } else {
	map {
	    $ModelPart_tabnames .= ", QModelParts d$_";
	} (1 .. $$study{LiabilityClassCnt});

        if (exists ($$study{TraitDist}) && $$study{TraitDist} eq 'T') {
            $colnames .= ", dXX.LittleSD, dXX.BigLittleSD";
            ($$study{ImprintingFlag} =~ /^y$/i) and $colnames .= ", dXX.LittleBigSD";
            $colnames .= ", dXX.BigSD";
        }
	$colnames .= ", dXX.LittleMean, dXX.BigLittleMean";
	($$study{ImprintingFlag} =~ /^y$/i) and $colnames .= ", dXX.LittleBigMean";
	$colnames .= ", dXX.BigMean";
    }
    ($$study{TraitType} eq 'QTT') and $colnames .= ", dXX.Threshold";
    
    map {
	$MPId_whereclause .= " and m.LC${_}MPId = d$_.MPId and m.LC${_}MPId = ?";
	($ModelPart_colnames .= $colnames) =~ s/XX/$_/g;
    } (1 .. $$study{LiabilityClassCnt});

    # We're guaranteed to have at least '[D|Q]ModelParts d1' in $ModelPart_tabnames, so it's
    # safe to refer to 'd1.DGF' in the static part of the SQL. Also, merlin_lk_prepare.pl
    # already excluded marker-only LK models (that is, DGF == -1), so we don't need to 
    # exclude those again here.

    ($sth = $dbh->prepare ("select m.ModelId, p.PedigreeSId, p.PedTraitPosCM, d1.DGF".
			   $ModelPart_colnames . " ".
			   "from Models m, PedigreePositions p".
			   $ModelPart_tabnames . " ".
			   "where p.StudyId = ? ".
			   "  and p.PedigreeSId in ".
			   "    (select PedigreeSId from ServerPedigrees ".
			   "     where ServerId = ?) ".
			   "  and p.PedPosId = m.PedPosId".
			   $MPId_whereclause))
	or die ("prepare select from Models failed, $DBI::errstr\n");
    print (ts(), "SQL is ". $sth->{Statement}. "\n");
    foreach (@MPIds) {
        print (ts(), "Selecting trait LK and combined LK models for MPIds ". join (', ', @$_). "\n");
	# Recall that each element in @MPIds is an arrayref containing one or more model part IDs
	$sth->execute ($$study{id}, $serverid, @$_)
	    or die ("execute select from Models failed, $DBI::errstr\n");
	@Models = @{$sth->fetchall_arrayref};
	(scalar (@Models) == 0) and die ("select returned no models\n");
	foreach $aref (@Models) {
	    ($modelid, $pedid, $traitpos) = splice (@$aref, 0, 3);
	    $key = join (':', @$aref);
	    if ($traitpos == -9999.99) {
		$$batch{$key}{$pedid}{trait} = $modelid;
	    } else {
		$$batch{$key}{$pedid}{$traitpos} = $modelid;
	    }
            push (@modelids, $modelid);
	}
	print (ts(), "Select returned ". scalar (@Models). " rows\n");
    }
    $sth->finish;
    db_update_modelids ($study, $serverid, \@modelids, 100);

    # This next bit only needs to be done by one process, so we're only going
    # to try it if 1) we're apparently running as a standalone (not an array) job,
    # or 2) if we are part of an array job, but we're running as task ID 1.

    if (defined ($ENV{SGE_TASK_ID}) && $ENV{SGE_TASK_ID} != 1) {
	print (ts(), "Skipping marker LK models\n");
	return ($batch);
    }

    print (ts(), "Selecting marker LK IDs\n");
    ($sth = $dbh->prepare ("select m.MarkerSetId, p.PedigreeSId, p.PedTraitPosCM ".
			   "from MarkerSetLikelihood m, PedigreePositions p ".
			   "where p.StudyId = ? ".
			   "  and p.PedigreeSId in ".
			   "    (select PedigreeSId from ServerPedigrees ".
			   "     where ServerId = ?) ".
			   "  and p.PedPosId = m.PedPosId ".
			   "  and m.ServerId is NULL ".
			   "  and m.Likelihood is NULL ".
			   "for update"))
	or die ("DBI prepare select marker LK IDs failed, $DBI::errstr\n");
    while (1) {
	if (! $sth->execute ($$study{id}, $serverid)) {
	    die ("DBI select marker LK models failed, $DBI::errstr\n");
	} elsif (defined ($DBI::errstr) && $DBI::errstr =~ /try restarting transaction/) {
	    print (ts(), "retrying select marker LK models\n");
	} else {
	    last;
	}
    }
    @Models = @{$sth->fetchall_arrayref};
    $sth->finish;
    print (ts(), "Select returned ". scalar (@Models). " rows\n");
    (scalar (@Models) == 0) and return ($batch);

    foreach $aref (@Models) {
        ($markersetid, $pedid, $traitpos) = @$aref;
	$$batch{marker}{$pedid}{$traitpos} = $markersetid;
        push (@markersetids, $markersetid);
    }
    db_update_markersetids ($study, $serverid, \@markersetids, 100);

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
    print (ts(), "Updating model part IDs in 1 batch of ", scalar (@$LGModelIds), "\n");
    while (1) {
	if ($sth->execute_array ({ArrayTupStatus => \@status}, $LGModelIds)) {
	    last;
	} elsif (! all_retry_errors (\@status, \@errors)) {
	    die ("execute update Models failed:\n", map { "$_\n" } @errors);
	}
	$retries++;
        ($debug) and print (ts(), "Updating model part IDs, total $retries retries\n");
    }
    ($auto || $dbh->commit) or die ("commit update LGModels failed, $DBI::errstr\n");
    $sth->finish;
    print (ts(), "Updated ". scalar (@$LGModelIds). " model part IDs with $retries retries\n");
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
    my $batchcount;

    db_connect ($study);
    ($sth = $dbh->prepare ("update Models set ServerId = $serverid, ".
			   "  MarkerCount = $$study{markercount}, ".
			   "  StartTime = NOW() ".
			   "where ModelId = ?"))
	or die ("DBI prepare update Models failed, $DBI::errstr\n");

    (defined ($batchsize)) or $batchsize = scalar (@$modelids);
    $batchcount = ceil (scalar (@$modelids) / $batchsize);
    print (ts(), "Updating trait LK and combined LK model ID in $batchcount batches of $batchsize records\n");
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
            ($debug) and print (ts(), "Updating trait LK and combined LK model IDs, on batch ", int ($count / $batchsize), ", total $retries retries\n");
	    $retries++;
	}
	($auto || $dbh->commit) or die ("commit update Models failed, $DBI::errstr\n");
    }
    $sth->finish;
    print (ts(), "Updated $count model IDs with $retries retries\n");
    return (1);
}


sub db_update_markersetids
{
    my ($study, $serverid, $markersetids, $batchsize) = @_;
    my $sth;
    my @status;
    my @batch;
    my @errors;
    my $count = 0;
    my $retries = 0;
    my $batchcount;

    db_connect ($study);
    ($sth = $dbh->prepare ("update MarkerSetLikelihood set ServerId = $serverid, ".
			   "  MarkerCount = $$study{markercount}, ".
			   "  StartTime = NOW() ".
			   "where MarkerSetId = ?"))
	or die ("DBI prepare update MarkerSetLikelihood failed, $DBI::errstr\n");

    (defined ($batchsize)) or $batchsize = scalar (@$markersetids);
    $batchcount = ceil (scalar (@$markersetids) / $batchsize);
    print (ts(), "Updating marker LK IDs in $batchcount batches of $batchsize records\n");
    while (scalar (@$markersetids)) {
        @batch = splice (@$markersetids, 0, $batchsize);
	while (1) {
	    if ($sth->execute_array ({ArrayTupleStatus => \@status}, \@batch)) {
		#print ("updated ", scalar (@batch), " MarkerSetLikelihoods\n");
		$count += scalar (@batch);
		last;
	    } elsif (! all_retry_errors (\@status, \@errors)) {
		die ("execute update MarkerSetLikelihood failed:\n", map { "$_\n" } @errors);
	    }
	    $retries++;
            ($debug) and print (ts(), "Updating marker LK IDs, on batch ", int ($count / $batchsize), ", total $retries retries\n");
	}
	($auto || $dbh->commit) or die ("commit update MarkerSetLikelihood failed, $DBI::errstr\n");
    }
    $sth->finish;
    print (ts(), "Updated $count marker LK IDs with $retries retries\n");
    return (1);
}


sub db_update_modelLKs
{
    my ($study, $modelids, $LKs, $batchsize) = @_;
    my $sth;
    my @status;
    my @idbatch;
    my @LKbatch;
    my @errors;
    my $count = 0;
    my $retries = 0;
    my $batchcount;

    db_connect ($study);
    ($sth = $dbh->prepare ("update Models ".
			   "set Likelihood = ?, EndTime = NOW() ".
			   "where ModelId = ?"))
	or die ("DBI prepare update of model LKs failed, $DBI::errstr\n");

    (defined ($batchsize)) or $batchsize = scalar (@$modelids);
    $batchcount = ceil (scalar (@$modelids) / $batchsize);
    print (ts(), "Updating model LKs in $batchcount batches of $batchsize records\n");
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
	    $retries++;
            ($debug) and print (ts(), "Updating model LKs, on batch ", int ($count / $batchsize), ", total $retries retries\n");
	}
	($auto || $dbh->commit) or die ("commit update model LKs failed, $DBI::errstr\n");
    }
    $sth->finish;
    print (ts(), "Updated $count model LKs with $retries retries\n");
    return (1);
}


sub db_update_markersetLKs
{
    my ($study, $markersetids, $LKs, $batchsize) = @_;
    my $sth;
    my @status;
    my @idbatch;
    my @LKbatch;
    my @errors;
    my $count = 0;
    my $retries = 0;
    my $batchcount;

    db_connect ($study);
    ($sth = $dbh->prepare ("update MarkerSetLikelihood ".
			   "set Likelihood = ?, EndTime = NOW() ".
			   "where MarkerSetId = ?"))
	or die ("DBI prepare update of model LKs failed, $DBI::errstr\n");

    (defined ($batchsize)) or $batchsize = scalar (@$markersetids);
    $batchcount = ceil (scalar (@$markersetids) / $batchsize);
    print (ts(), "Updating marker LKs in $batchcount batches of $batchsize records\n");
    while (scalar (@$markersetids)) {
	@idbatch = splice (@$markersetids, 0, $batchsize);
	@LKbatch = splice (@$LKs, 0, $batchsize);
	while (1) {
	    if ($sth->execute_array ({ArrayTupleStatus => \@status}, \@LKbatch, \@idbatch)) {
		#print ("updated ", scalar (@idbatch), " markerset LKs\n");
		$count += scalar (@idbatch);
		last;
	    } elsif (! all_retry_errors (\@status, \@errors)) {
		die ("execute update markerset LKs failed:\n", map { "$_\n" } @errors);
	    }
	    $retries++;
            ($debug) and print (ts(), "Updating marker LKs, on batch ", int ($count / $batchsize), ", total $retries retries\n");
	}
	($auto || $dbh->commit) or die ("commit update markerset LKs failed, $DBI::errstr\n");
    }
    $sth->finish;
    print (ts(), "Updated $count marker LKs with $retries retries\n");
    return (1);
}


sub db_set_server_status
{
    my ($study, $serverid, $status) = @_;
    my $sth;

    db_connect ($study);
    (($sth = $dbh->prepare ("call ServerSignOff (?, ?)")) &&
     $sth->execute ($serverid, $status))
	or warn ("DBI call ServerSignOff() failed, $DBI::errstr\n");
    ($auto || $dbh->commit) or warn ("DBI commit update Servers failed, $DBI::errstr\n");
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
