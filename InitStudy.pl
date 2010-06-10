#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use DBI; # Database interaction
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
KelvinDataset->VERSION (1.00);
require KelvinConfig;
KelvinConfig->VERSION (1.00);
require KelvinFamily;
KelvinFamily->VERSION (1.00);

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

sub perform_study
{
    my ($config) = @_;

    $ {$config->isConfigured ("Study")}[0] =~ /(\d+)\s+(\w+)\s+\"(.+)\"\s+(\w+)\s+(\w+)/;
    my $StudyId = $1; my $StudyRole = lc($2); my $DBIDatabase = $3; my $Username = $4; my $Password = $5;

    my $dbh;
    until ($dbh = DBI->connect("dbi:$DBIDatabase", $Username, $Password,
			       { RaiseError => 0, PrintError => 0, AutoCommit => 1 })) {
	sleep(5);
	warn "Cannot connect to $DBIDatabase: $DBI::errstr, retrying!";
    }

    my $ReferenceMapId; my $LiabilityClassCnt; my $ImprintingFlag;

    # We always want a full map

    my $href = {};
    my $dataset;

    # Parse in the data files...
    $$href{PedigreeFile} = $ {$config->isConfigured ("PedigreeFile")}[0];
    $$href{LocusFile} = $ {$config->isConfigured ("LocusFile")}[0];
    $$href{MapFile} = $ {$config->isConfigured ("MapFile")}[0];
    $dataset = KelvinDataset->new ($href)
	or error ("KelvinDataset->new failed");
    $ReferenceMapId = find_or_insert_map($$href{MapFile}, $dataset, $dbh);

    # Get the Studies row indicated by the configuration...
    my $sth = $dbh->prepare("Select ReferenceMapId, LiabilityClassCnt, ImprintingFlag from Studies where StudyId = ?")
	or die "Couldn't prepare Studies selection: $dbh->errstr";
    $sth->execute($StudyId) or die "Couldn't execute Studies selection: $dbh->errstr";
    my @Results;
    if (@Results = $sth->fetchrow_array()) {
	# Check map, liability classes and imprinting
	my $StudyMapId = $Results[0]; $LiabilityClassCnt = $Results[1]; $ImprintingFlag = lc($Results[2]);
	error ("Mismatch in map")
	    if ($StudyMapId ne $ReferenceMapId);
	if ($config->isConfigured ("LiabilityClasses")) {
	    error ("Mismatch in liability study/run liability classes")
		if ($ {$config->isConfigured ("LiabilityClasses")}[0] ne $LiabilityClassCnt);
	} else {
	    error ("LiabilityClasses directive missing and Studies row specifies $LiabilityClassCnt classes")
		if ($LiabilityClassCnt ne 1);
	}
	if ($config->isConfigured ("Imprinting")) {
	    error ("Mismatch in imprinting flag ('on' in directives, 'off' in Studies row)") if ($ImprintingFlag ne "y");
	} else {
	    error ("Mismatch in imprinting flag ('off' in directives, 'on' in Studies row)") if ($ImprintingFlag eq "y");
	}
	$sth->finish;
    } else {
	error ("Client run must setup Studies and other tables before servers can run")
	    if ($StudyRole ne "client");
	$LiabilityClassCnt = $ {$config->isConfigured ("LiabilityClasses")}[0];
	$ImprintingFlag = $config->isConfigured ("Imprinting") ? 'Y' : 'N';

	# Insert Studies row with information from this configuration
	$dbh->do("Insert into Studies (StudyId, ReferenceMapId, LiabilityClassCnt, ImprintingFlag) values (?,?,?,?)",
		 undef, $StudyId, $ReferenceMapId, $LiabilityClassCnt, $ImprintingFlag)
	    or die "Cannot insert Studies row: $DBI::errstr";
    }

    # Make sure all PedigreeSIds and SingleModelTimes are present

    # Make a list of PedigreeSIds by reading the pedigree file
    while (my $ped = $dataset->readFamily) { 
#	print Dumper($ped);
	my $PedigreeSId = $$ped{pedid};
	# Don't care if this fails...
	$dbh->do("Insert into Pedigrees (StudyId, PedigreeSId) values (?,?)",
		 undef, $StudyId, $PedigreeSId);
    }
    # 'Freshen' the SingleModelTimes table by creating a a replacement...
    $dbh->do("Create temporary table SMTs Select a.StudyId, b.PedigreeSId, c.ChromosomeNo, d.MarkerName, 2 ".
	     "from Studies a, Pedigrees b, Markers c, MapMarkers d where ".
	     "a.StudyId = 17 AND b.StudyId = a.StudyId AND c.Name = d.MarkerName AND d.MapId = a.ReferenceMapId");
    # ...and then inserting any new combinations.
    $dbh->do("Insert into SingleModelTimes (StudyId, PedigreeSId, ChromosomeNo, MarkerName) ".
	     "Select a.StudyId, a.PedigreeSId, a.ChromosomeNo, a.MarkerName from ".
	     "SMTs a left outer join SingleModelTimes b on ".
	     "a.StudyId = b.StudyId AND a.PedigreeSId = b.PedigreeSId AND a.MarkerName = b.MarkerName ".
	     "where b.StudyId IS NULL");

    # 'Freshen' the Positions table
    # Three cases we know of: marker, individual value, and range specification, and all can be in lists
    my $JointTPs = join(',', @{$config->isConfigured ("TraitPositions")});
    $JointTPs =~ s/\s+//g;
    my @TPs = split(',',$JointTPs);
    my $ChromosomeNo = $$dataset{chromosome};
    for my $TP (@TPs) {
	if ($TP eq "marker") {
	    $dbh->do("Insert into Positions (StudyId, ChromosomeNo, RefTraitPosCM) ".
		     "Select $StudyId, $ChromosomeNo, AveragePosCM from ".
		     "MapMarkers where MapId = $ReferenceMapId AND AveragePosCM NOT in ".
		     "(Select distinct RefTraitPosCM from Positions where ".
		     "StudyId = $StudyId)", undef);
	} elsif ($TP =~ /(\d*.?\d*)-(\d*.?\d*):(\d*.?\d*)/) {
	    my $PosCM = $1;
	    do {
		$dbh->do("Insert into Positions (StudyId, ChromosomeNo, RefTraitPosCM) values (?,?,?)",
			 undef, $StudyId, $ChromosomeNo, $PosCM);
		$PosCM += $2;
	    } while ($PosCM <= $3);
	} else {
	    $dbh->do("Insert into Positions (StudyId, ChromosomeNo, RefTraitPosCM) values (?,?,?)",
		     undef, $StudyId, $ChromosomeNo, $TP);
	}
    }

    # Finally, freshen the PedigreePositions as needed...
    my $MPMarkers = $ {$config->isConfigured ("Multipoint")}[0];
    print "Running with $MPMarkers multipoint markers\n";
    
    $dbh->do("Create temporary table PPs Select a.StudyId, a.PedigreeSId, b.ChromosomeNo, ".
	     "b.RefTraitPosCM, $MPMarkers 'MarkerCount' ".
	     "from Pedigrees a, Positions b where a.StudyId = b.StudyId");
    $dbh->do("Insert into PedigreePositions (StudyId, PedigreeSId, ChromosomeNo, RefTraitPosCM, MarkerCount) ".
	     "Select a.StudyId, a.PedigreeSId, a.ChromosomeNo, a.RefTraitPosCM, a.MarkerCount from ".
	     "PPs a left outer join PedigreePositions b on ".
	     "a.StudyId = b.StudyId AND a.PedigreeSId = b.PedigreeSId AND a.ChromosomeNo = b.ChromosomeNo AND ".
	     "a.RefTraitPosCM = b.RefTraitPosCM AND a.MarkerCount = b.MarkerCount ".
	     "where b.StudyId IS NULL");
}

sub find_or_insert_map
{
    my $MapFile = shift();
    my $dataset = shift();
    my $dbh = shift();

    my $MapScale = uc(substr $$dataset{mapfunction},0,1);
    my $MapId = 0;
	
#    print "fields are ".Dumper($$dataset{mapfields})."\n";

    # Get the Maps row...
    my $sth = $dbh->prepare("Select MapId from Maps where MapScale = ? AND Description like ?")
	or die "Couldn't prepare Maps selection: $dbh->errstr";
    $sth->execute($MapScale, $MapFile."%") or die "Couldn't execute Maps selection: $dbh->errstr";
    my @Results;
    if (@Results = $sth->fetchrow_array()) {
	# Got it, return it...
	$MapId = $Results[0];
	$sth->finish;
    } else {
	# No such map, insert it
	$dbh->do("Insert into Maps (MapScale, Description) values (?,?)", undef, $MapScale, $MapFile)
	    or die "Cannot insert Maps row: $DBI::errstr";
	$sth->finish;
	$sth = $dbh->prepare("Select LAST_INSERT_ID()");
	$sth->execute;
	$sth->bind_col(1, \$MapId);
	$sth->fetch;
	$sth->finish;
    }
    # Ensure that all markers are present even if they're a superset of an existing map
    my %markers = %{$$dataset{markers}};
    my @markerOrder = @{$$dataset{markerorder}};
    my $ChromosomeNo = $$dataset{chromosome};
    # Find or insert the markers. Marker names must be identical between maps for interpolation.
    foreach (@markerOrder) {
	# Don't care if this fails...
	$dbh->do("Insert into Markers (Name, ChromosomeNo) values (?,?)", undef, $_, $ChromosomeNo);
	$dbh->do("Insert into MapMarkers (MapId, MarkerName, AveragePosCM) values (?,?,?)",
		 undef, $MapId, $_, $markers{$_}{avgpos});
    }
    return $MapId;
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

