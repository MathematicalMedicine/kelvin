#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use DBI; # Database interaction
use DBI qw(:sql_types);
$|=1; # Show the output when I say so.

my $KELVIN_ROOT='NO_KELVIN_ROOT';

my $usage = "usage: $0 <configfile> [--directive ... ]\n";
my $config;
my $configFile;
my ($directive, $args);
my $arg;
my $idx = 0;
my $debug = 0;

# For all these, we allow environment variables to override everything, even if
# values were set during installation.
if ($ENV{KELVIN_ROOT}) {
    ($KELVIN_ROOT !~ /no_kelvin_root/i)
	and warner ("overriding installed KELVIN_ROOT with '$ENV{KELVIN_ROOT}' from environment");
    $KELVIN_ROOT = $ENV{KELVIN_ROOT};
} elsif ($KELVIN_ROOT =~ /no_kelvin_root/i) {
    $KELVIN_ROOT = dirname ($0);
    warner ("no KELVIN_ROOT defined by installation, using '$KELVIN_ROOT'");
}

unshift (@INC, $KELVIN_ROOT);
require KelvinDataset;
KelvinDataset->VERSION (1.40);
require KelvinConfig;
KelvinConfig->VERSION (1.20);
require KelvinFamily;
KelvinFamily->VERSION (1.40);

($configFile = shift (@ARGV))
    or die ($usage);
($config = KelvinConfig->new ($configFile))
    or error ("new KelvinConfig failed: $KelvinConfig::errstr");
if (@ARGV) {
    while (defined ($directive = $ARGV[$idx++])) {
	($directive =~ s/^--//) or die ($usage);
	$args = undef;
	while (defined ($ARGV[$idx]) && $ARGV[$idx] !~ /^--/) {
	    $arg = $ARGV[$idx++];
	    $args = (defined ($args)) ? $args .= " $arg" : $arg;
	}
	($config->addDirective ($directive, $args))
	    or error ("$KelvinConfig::errstr on command line");
    }
    $config->validate
	or error ("$KelvinConfig::errstr");
}
($debug) and print Dumper ($config);

($config->isConfigured ("Study"))
    or error ("Study directive must be provided");
($config->isConfigured ("Multipoint"))
    or error ("Must be a multipoint analysis");

perform_study($config);

#
# Current limitations:
#
# 1. There can be only one genotype map per pedigree. Otherwise we have to track which position
#    ranges to use for each map for each pedigree and it gets very confusing.
# 2. Multiple disjoint client runs will produce *cross-products* of requested positions and 
#    pedigrees whether you want that or not, i.e. client A requesting ped 1 & 2 for pos 1-10
#    followed by a run of client B requesting ped 2 and 3 for pos 11-20 will also end up requesting
#    pos 11-20 for ped 1 and pos 1-10 for ped 3. To avoid this I would have to associate both
#    peds and pos with some kind of generated client ID to keep them together during the populating
#    joins. Not horribly difficult, but tough enough to not warrant it when it is unlikely to be
#    needed.
#    

sub perform_study
{
    my ($config) = @_;

    # The following regex allows quoting of StudyLabel and Password, which bears narration: optionally
    # match a single or full quote, then a string excluding single or full quotes, and finally, if
    # some quote was found, require its closure. Only shortcoming is not allowing embedded single or
    # full quotes.
    ${$config->isConfigured ("Study")}[0] =~ /(["'])?([^"']+)(?(1)\1|)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(["'])?([^"']+)(?(7)\7|)\s+(\S+)\s+(\S+)/;

    my $StudyLabel = $2; my $StudyRole = uc($3); my $DBIHost = $4; my $DBIDatabase = $5; my $Username = $6; my $Password = $8; my $PedigreeRegEx = $9; my $PedigreeNotRegEx = $10;
    my $DBIConnectionString = "mysql:host=".$DBIHost.":database=$DBIDatabase";
    my $MapId; my $LiabilityClasses = 1; my $ImprintingFlag = 'n';

    $LiabilityClasses = $ {$config->isConfigured ("LiabilityClasses")}[0]
	if ($config->isConfigured ("LiabilityClasses"));
    $ImprintingFlag = 'y' if ($config->isConfigured ("Imprinting"));

    my $dbh;
    until ($dbh = DBI->connect("dbi:$DBIConnectionString", $Username, $Password,
			       { RaiseError => 1, PrintError => 1, AutoCommit => 1 })) {
	sleep(5);
	warn "Cannot connect to $DBIDatabase: $DBI::errstr, retrying!";
    }

    # Parse in the data files...

    my $href = {};
    my $dataset;
    $$href{PedigreeFile} = $ {$config->isConfigured ("PedigreeFile")}[0];
    $$href{LocusFile} = $ {$config->isConfigured ("LocusFile")}[0];
    $$href{MapFile} = $ {$config->isConfigured ("MapFile")}[0];
    $dataset = KelvinDataset->new ($href)
	or error ("KelvinDataset->new failed, $KelvinDataset::errstr");
 
    # Set the undefined phenocode so KelvinFamily can identify unphenotyped individuals
    my $aref = $config->isConfigured ("PhenoCodes");
    $dataset->setUndefPhenocode ((split (/,\s*/, $$aref[0]))[0]);

    # Get or insert the Study
    my $StudyId;
    #my $sth = $dbh->prepare("Call GetStudyId(?, ?, ?, ?)")
    #	or die "Couldn't prepare GetStuydId: $dbh->errstr";
    # $sth->bind_param_inout(4, \$StudyId, 11, {TYPE=>SQL_INTEGER});
    # $sth->execute($StudyName, $LiabilityClasses, $ImprintingFlag) or die "Couldn't execute GetStudyId: $dbh->errstr";
    # The above nets me: DBD::mysql::st bind_param_inout failed: Output parameters not implemented at ../InitStudy.pl line 121,
    # so the workaround is:
    $dbh->do("Call GetStudyId(\'$StudyLabel\', $LiabilityClasses, \'$ImprintingFlag\', \@StudyId)");
    $StudyId = $dbh->selectrow_array('Select @StudyId');

    # Get or insert the Map
    my $MapScale = uc(substr $$dataset{mapfunction},0,1);
    my $MapFile = $$href{MapFile};
    my $MapDescription;
    $dbh->do("Call GetMapId(\'$StudyId\', \'".(($StudyRole eq "CLIENT")? 1:0)."\', \'$MapScale\', \'$MapFile\', \@MapId, \@MapScale, \@MapDescription)");
    ($MapId, $MapScale, $MapDescription) = $dbh->selectrow_array('Select @MapId, @MapScale, @MapDescription');

    # Ensure that all markers are present even if they're a superset of an existing map
    my @markerOrder = @{$$dataset{markerorder}};
    my %markers = %{$$dataset{markers}};
    my $ChromosomeNo = $$dataset{chromosome};
    # Find or insert the markers. Marker names must be identical between maps for interpolation.
    foreach (@markerOrder) {
	# Don't care if this fails...
	$dbh->do("Insert ignore into Markers (StudyId, Name, ChromosomeNo) values (?,?,?)", undef, $StudyId, $_, $ChromosomeNo);
	$dbh->do("Insert ignore into MapMarkers (StudyId, MapId, MarkerName, AvePosCM) values (?,?,?,?)",
		 undef, $StudyId, $MapId, $_, $markers{$_}{avgpos});
    }

    # Make sure all PedigreeSIds and SingleModelTimes are present

    # Get list of PedigreeSIds by reading the pedigree file. Add if we're a client, update if we're a server

    my $ped;
    while ($ped = $dataset->readFamily) { 
#	print Dumper($ped);
	my $PedigreeSId = $$ped{pedid};
	if ($StudyRole eq "CLIENT") {
	    # Client -- just slam 'em in, don't care if this fails with duplicates...NOTE that
            # Inclusion and exclusion regexps are ignored since the kelvin client processes
            # all pedigrees. Only servers have a restricted view.
	    $dbh->do("Insert ignore into Pedigrees (StudyId, PedigreeSId) values (?,?)",
		     undef, $StudyId, $PedigreeSId);
	} else {
	    # Server -- no worries about errors here either...
	    
	    $dbh->do("Update ignore Pedigrees set GenotypeMapId = ? where StudyId = ? AND PedigreeSId = ?",
		     undef, $MapId, $StudyId, $PedigreeSId)
		if (($PedigreeSId =~ $PedigreeRegEx) && ($PedigreeSId !~ $PedigreeNotRegEx));
	}
    }
    (! defined ($ped)) and error ($KelvinDataset::errstr);

    if ($StudyRole eq "SERVER") {
	$dbh->do("call BadScaling(?)", undef, $StudyId);
	# Add the following to hold trait likelihood references
	$dbh->do("Insert ignore into PedigreePositions ".
		 "(StudyId, PedigreeSId, ChromosomeNo, RefTraitPosCM, PedTraitPosCM, MarkerCount, FreeModels) ".
		 "Select distinct $StudyId, PedigreeSId, ChromosomeNo, -9999.99, -9999.99, MarkerCount, 0 from PedigreePositions where StudyId = $StudyId", undef);

	return; # All client-only from here on out...
    }

    # 'Freshen' the Positions tables
    # Three cases we know of: marker, individual value, and range specification, and all can be in lists
    my $JointTPs = join(',', @{$config->isConfigured ("TraitPositions")});
    $JointTPs =~ s/\s+//g;
    my @TPs = split(',',$JointTPs);
    $ChromosomeNo = $$dataset{chromosome};
    # First find where "begin" and "end" are for the current set of markers.
    my $lastMarkerPos = 0;
    my $firstMarkerPos = 9999.9;
    %markers = %{$$dataset{markers}};
    foreach (keys %markers) {
	$lastMarkerPos = $markers{$_}{avgpos} if ($markers{$_}{avgpos} > $lastMarkerPos);
	$firstMarkerPos = $markers{$_}{avgpos} if ($markers{$_}{avgpos} < $firstMarkerPos);
    }
    for my $TP (@TPs) {
	if ($TP =~ /.*-end:(\d*.?\d*)/) {
	    my $newEnd = $lastMarkerPos + $1;
	    $TP =~ s/-end:/-$newEnd:/;
	}
	if ($TP =~ /begin-(\d*.?\d*):(\d*.?\d*)/) {
	    my $newBegin = int($firstMarkerPos / $2) * $2;
	    $TP =~ s/begin-/$newBegin-/;
	}
	if (uc $TP eq "MARKER") {
	    $dbh->do("Insert ignore into Positions (StudyId, ChromosomeNo, RefTraitPosCM) ".
		     "Select $StudyId, $ChromosomeNo, AvePosCM from ".
		     "MapMarkers where StudyId = $StudyId AND MapId = $MapId AND AvePosCM NOT in ".
		     "(Select distinct RefTraitPosCM from Positions where ".
		     "StudyId = $StudyId)", undef);
	} elsif ($TP =~ /(\d*.?\d*)-(\d*.?\d*):(\d*.?\d*)/) {
	    my ($begin, $end, $inc, $precision, $va) = ($1, $2, $3, $3, 0);
	    $precision =~ s/^[^\.]\.?(\d*)?/$1/;
	    my $PosCM;
	    do {
		$PosCM = sprintf ("%0.*f", length ($precision), $begin + $inc * $va++);
		$dbh->do("Insert ignore into Positions (StudyId, ChromosomeNo, RefTraitPosCM) values (?,?,?)",
			 undef, $StudyId, $ChromosomeNo, $PosCM);
	    } while ($PosCM <= $end);
	} else {
	    $dbh->do("Insert ignore into Positions (StudyId, ChromosomeNo, RefTraitPosCM) values (?,?,?)",
		     undef, $StudyId, $ChromosomeNo, $TP);
	}
    }

    # Finally, freshen the PedigreePositions as needed...
    my $MPMarkers = $ {$config->isConfigured ("Multipoint")}[0];
    
    $dbh->do("Create temporary table PPs Select a.StudyId, a.PedigreeSId, b.ChromosomeNo, ".
	     "b.RefTraitPosCM, $MPMarkers 'MarkerCount' ".
	     "from Pedigrees a, Positions b where a.StudyId = $StudyId AND a.StudyId = b.StudyId");
    $dbh->do("Insert into PedigreePositions (StudyId, PedigreeSId, ChromosomeNo, RefTraitPosCM, MarkerCount) ".
	     "Select a.StudyId, a.PedigreeSId, a.ChromosomeNo, a.RefTraitPosCM, a.MarkerCount from ".
	     "PPs a left outer join PedigreePositions b on ".
	     "a.StudyId = $StudyId AND a.StudyId = b.StudyId AND ".
	     "a.PedigreeSId = b.PedigreeSId AND ".
	     "a.ChromosomeNo = b.ChromosomeNo AND ".
	     "a.RefTraitPosCM = b.RefTraitPosCM ".
	     "where b.StudyId IS NULL");
    $dbh->do("call BadScaling(?)", undef, $StudyId);

    return;
}

sub fatal
{
    die ("FATAL - ABORTING, @_");
}

sub error
{
    die ("ERROR - EXITING, @_\n");
}

sub warner
{
    warn ("WARNING, @_\n");
}

