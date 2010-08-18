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

    $ {$config->isConfigured ("Study")}[0] =~ /(\d+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\S+)/;
    my $StudyId = $1; my $StudyRole = lc($2); my $DBIHost = $3; my $DBIDatabase = $4; my $Username = $5; my $Password = $6; my $PedigreeRegEx = $7;
    my $DBIConnectionString = "mysql:host=$DBIHost:database=$DBIDatabase";
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

    # We always want a full map whether it is for the client (ReferenceMap) or server (GenotypeMaps)

    $MapId = find_or_insert_map($$href{MapFile}, $dataset, $dbh);

    # Get the Studies row indicated by the configuration...

    my $sth = $dbh->prepare("Select ReferenceMapId, LiabilityClassCnt, ImprintingFlag from Studies where StudyId = ?")
	or die "Couldn't prepare Studies selection: $dbh->errstr";
    $sth->execute($StudyId) or die "Couldn't execute Studies selection: $dbh->errstr";

    if (my @Results = $sth->fetchrow_array()) {

	# Study exists, check map (if client), liability classes and imprinting

	error ("Mismatch in map")
	    if ($StudyRole eq "client" and $Results[0] ne $MapId);
	error ("Mismatch in liability classes")
	    if ($Results[1] ne $LiabilityClasses);
	error ("Mismatch in imprinting flag")
	    if (lc($Results[2]) ne $ImprintingFlag);
	$sth->finish;

    } else {

	# Study does not exist, OK for the client, bad for a server.

	error ("Client run must setup Studies and other tables before servers can run")
	    if ($StudyRole ne "client");

	# Insert Studies row with information from this configuration
	$dbh->do("Insert into Studies (StudyId, ReferenceMapId, LiabilityClassCnt, ImprintingFlag) values (?,?,?,?)",
		 undef, $StudyId, $MapId, $LiabilityClasses, $ImprintingFlag)
	    or die "Cannot insert Studies row: $DBI::errstr";
    }

    # Make sure all PedigreeSIds and SingleModelTimes are present

    # Get list of PedigreeSIds by reading the pedigree file. Add if we're a client, update if we're a server

    my $ped;
    while ($ped = $dataset->readFamily) { 
#	print Dumper($ped);
	my $PedigreeSId = $$ped{pedid};
	if ($StudyRole eq "client") {
	    # Client -- just slam 'em in, don't care if this fails with duplicates...
	    $dbh->do("Insert into Pedigrees (StudyId, PedigreeSId) values (?,?)",
		     undef, $StudyId, $PedigreeSId);
	} else {
	    # Server -- no worries about errors here either...
	    $dbh->do("Update Pedigrees set GenotypeMapId = ? where PedigreeSId = ?",
		     undef, $MapId, $PedigreeSId);
	}
    }
    (! defined ($ped)) and error ($KelvinDataset::errstr);

    if ($StudyRole eq "server") {
	$dbh->do("call BadScaling(?)", undef, $StudyId);
	return; # All client-only from here on out...
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

    # 'Freshen' the Positions tables
    # Three cases we know of: marker, individual value, and range specification, and all can be in lists
    my $JointTPs = join(',', @{$config->isConfigured ("TraitPositions")});
    $JointTPs =~ s/\s+//g;
    my @TPs = split(',',$JointTPs);
    my $ChromosomeNo = $$dataset{chromosome};
    for my $TP (@TPs) {
	if ($TP eq "marker") {
	    $dbh->do("Insert into Positions (StudyId, ChromosomeNo, RefTraitPosCM) ".
		     "Select $StudyId, $ChromosomeNo, AvePosCM from ".
		     "MapMarkers where MapId = $MapId AND AvePosCM NOT in ".
		     "(Select distinct RefTraitPosCM from Positions where ".
		     "StudyId = $StudyId)", undef);
	} elsif ($TP =~ /(\d*.?\d*)-(\d*.?\d*):(\d*.?\d*)/) {
	    my $PosCM = $1;
	    do {
		$dbh->do("Insert into Positions (StudyId, ChromosomeNo, RefTraitPosCM) values (?,?,?)",
			 undef, $StudyId, $ChromosomeNo, $PosCM);
#		$PosCM = sprintf ("%.2f", $PosCM + $3);
		$PosCM += $3;
	    } while ($PosCM <= $2);
	} else {
	    $dbh->do("Insert into Positions (StudyId, ChromosomeNo, RefTraitPosCM) values (?,?,?)",
		     undef, $StudyId, $ChromosomeNo, $TP);
	}
    }

    # Finally, freshen the PedigreePositions as needed...
    my $MPMarkers = $ {$config->isConfigured ("Multipoint")}[0];
    
    $dbh->do("Create temporary table PPs Select a.StudyId, a.PedigreeSId, b.ChromosomeNo, ".
	     "b.RefTraitPosCM, $MPMarkers 'MarkerCount' ".
	     "from Pedigrees a, Positions b where a.StudyId = b.StudyId");
    $dbh->do("Insert into PedigreePositions (StudyId, PedigreeSId, ChromosomeNo, RefTraitPosCM, MarkerCount) ".
	     "Select a.StudyId, a.PedigreeSId, a.ChromosomeNo, a.RefTraitPosCM, a.MarkerCount from ".
	     "PPs a left outer join PedigreePositions b on ".
	     "a.StudyId = b.StudyId AND a.PedigreeSId = b.PedigreeSId AND a.ChromosomeNo = b.ChromosomeNo AND ".
	     "a.RefTraitPosCM = b.RefTraitPosCM ".
	     "where b.StudyId IS NULL");
    $dbh->do("call BadScaling(?)", undef, $StudyId);

    return;
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
    if (my @Results = $sth->fetchrow_array()) {
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
    $dbh->{PrintError} = 0; $dbh->{RaiseError} = 0;
    foreach (@markerOrder) {
	# Don't care if this fails...
	$dbh->do("Insert into Markers (Name, ChromosomeNo) values (?,?)", undef, $_, $ChromosomeNo);
	$dbh->do("Insert into MapMarkers (MapId, MarkerName, AvePosCM) values (?,?,?)",
		 undef, $MapId, $_, $markers{$_}{avgpos});
    }
    $dbh->{PrintError} = 1; $dbh->{RaiseError} = 1;
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

