#!/usr/bin/env perl

# Creates a new Study within a LKS database and adds that which the study will
# be analyzing to same.

package LKS_InitStudy;

our $VERSION = '$Id$';

use strict;
use warnings FATAL => qw(all);
use lib ($ENV{TOOLPATH} ? split (/:/, $ENV{TOOLPATH}) : ' ');
use English qw( -no_match_vars );

use DBI; # Database interaction
use DBI qw(:sql_types);

use RunSQLScript;
use KelvinDataset 1.40;
use KelvinConfig 1.20;
use KelvinFamily 1.40;

use Exporter qw(import);
our @EXPORT_OK = qw(perform_study);

$OUTPUT_AUTOFLUSH= 1 ;  # Show the output when I say so.

my $optssetup = [
    "Creates a new Study within a LKS database and adds that which the study "
            . "will be analyzing to same.",
    "conffile!",
];

sub main {
    $main::VERSION = $VERSION;  # needed for --version option to work
    init_study(\@ARGV);
}

sub init_study {
    my $opts = BCMM::CLIParser::autoparse_cmdline($optssetup, \@ARG);
    
    my $config = KelvinConfig->new($$opts{conffile})
            or error("new KelvinConfig failed: $KelvinConfig::errstr");
    
    $config->validate or error("$KelvinConfig::errstr");
    $config->isConfigured("Study")
            or error("Study directive must be provided");
    $config->isConfigured("Multipoint")
            or error("Must be a multipoint analysis");
    
    perform_study($config);
}

#
# Current limitations:
#
# 1. There can be only one genotype map per pedigree. Otherwise we have to
#    track which position ranges to use for each map for each pedigree and it
#    gets very confusing.
# 2. Multiple disjoint client runs will produce *cross-products* of requested
#    positions and pedigrees whether you want that or not, i.e. client A
#    requesting ped 1 & 2 for pos 1-10 followed by a run of client B requesting
#    ped 2 and 3 for pos 11-20 will also end up requesting pos 11-20 for ped 1
#    and pos 1-10 for ped 3. To avoid this I would have to associate both
#    peds and pos with some kind of generated client ID to keep them together
#    during the populating joins. Not horribly difficult, but tough enough to
#    not warrant it when it is unlikely to be needed.
#    

sub perform_study {
    my ($config) = @ARG;
    
    my $study = $config->readStudyLine();
    
    error("StudyLabel ($$study{label}) in STUDY directive too long")
            if (length($$study{label}) > 63);
    my $MapId; my $LiabilityClasses = 1; my $ImprintingFlag = 'n';
    
    $LiabilityClasses = $ {$config->isConfigured ("LiabilityClasses")}[0]
            if ($config->isConfigured ("LiabilityClasses"));
    $ImprintingFlag = 'y' if ($config->isConfigured ("Imprinting"));
    
    my $dbh = RunSQLScript::connect_to_db($$study{host}, $$study{database},
            $$study{username}, $$study{password});
    $dbh->{PrintError} = 1;
            # we also get RaiseError and AutoCommit for free from
            # connect_to_db()
    
    # Parse in the data files...
    
    my $datafiles = {};
    my $dataset;
    $$datafiles{PedigreeFile} = ${$config->isConfigured("PedigreeFile")}[0];
    $$datafiles{LocusFile} = ${$config->isConfigured("LocusFile")}[0];
    $$datafiles{MapFile} = ${$config->isConfigured("MapFile")}[0];
    $dataset = KelvinDataset->new($datafiles)
            or error("KelvinDataset->new failed, $KelvinDataset::errstr");
 
    # Set the undefined phenocode so KelvinFamily can identify unphenotyped
    # individuals
    my $phenocodes = $config->isConfigured("PhenoCodes");
    $dataset->setUndefPhenocode((split (/,\s*/, $$phenocodes[0]))[0]);
    
    # Get or insert the Study
    my $StudyId;
    #my $sth = $dbh->prepare("Call GetStudyId(?, ?, ?, ?)")
    #        or die "Couldn't prepare GetStuydId: $dbh->errstr";
    #$sth->bind_param_inout(4, \$StudyId, 11, {TYPE=>SQL_INTEGER});
    #$sth->execute($StudyName, $LiabilityClasses, $ImprintingFlag)
    #        or die "Couldn't execute GetStudyId: $dbh->errstr";
    # The above nets me:
    # DBD::mysql::st bind_param_inout failed: Output parameters not implemented
    # at ../InitStudy.pl line 121
    # so the workaround is:
    $dbh->do("CALL GetStudyId(\'$$study{label}\', $LiabilityClasses, "
            . "\'$ImprintingFlag\', \@StudyId)");
    $StudyId = $dbh->selectrow_array('SELECT @StudyId');
    
    # Get or insert the Map
    my $MapScale = uc(substr $$dataset{mapfunction}, 0, 1);
    my $MapFile = $$datafiles{MapFile};
    my $MapDescription;
    $dbh->do("CALL GetMapId(\'$StudyId\', \'" . (($$study{role} eq "CLIENT")? 1:0)
            . "\', \'$MapScale\', \'$MapFile\', \@MapId, \@MapScale, "
            . "\@MapDescription)");
    ($MapId, $MapScale, $MapDescription) = $dbh->selectrow_array(
            'SELECT @MapId, @MapScale, @MapDescription');
    
    # Ensure that all markers are present even if they're a superset of an
    # existing map
    my @markerOrder = @{$$dataset{markerorder}};
    my %markers = %{$$dataset{markers}};
    my $ChromosomeNo = $$dataset{chromosome};
    # Find or insert the markers. Marker names must be identical between maps
    # for interpolation.
    foreach my $marker (@markerOrder) {
        # Don't care if this fails...
        $dbh->do("INSERT IGNORE INTO Markers (StudyId, Name, ChromosomeNo) "
                . "VALUES (?,?,?)", undef, $StudyId, $marker, $ChromosomeNo);
        $dbh->do("INSERT IGNORE INTO MapMarkers (StudyId, MapId, MarkerName, "
                . "AvePosCM) VALUES (?,?,?,?)",
                undef, $StudyId, $MapId, $marker, $markers{$marker}{avgpos});
    }
    
    # Make sure all PedigreeSIds and SingleModelTimes are present
    
    # Get list of PedigreeSIds by reading the pedigree file. Add if we're a
    # client, update if we're a server
    
    my $ped;
    while ($ped = $dataset->readFamily) { 
        #print Dumper($ped);
        my $PedigreeSId = $$ped{pedid};
        if ($$study{role} eq "CLIENT") {
            # Client -- just slam 'em in, don't care if this fails with
            # duplicates...
            # NOTE that inclusion and exclusion regexps are ignored since the
            # kelvin client processes all pedigrees. Only servers have a
            # restricted view.
            # FIXME: possibly should explicitly check errors and verify that
            # the only errors ARE duplicates, altho INSERT IGNORE *should* take
            # care of that for us...
            $dbh->do("INSERT IGNORE INTO Pedigrees (StudyId, PedigreeSId) "
                    . "VALUES (?,?)", undef, $StudyId, $PedigreeSId);
        } else {
            # Server -- no worries about duplication errors here either...
            if (($PedigreeSId =~ $$study{pedregex})
                    && ($PedigreeSId !~ $$study{pednotregex})) {
                $dbh->do("UPDATE IGNORE Pedigrees SET GenotypeMapId = ? WHERE "
                        . "StudyId = ? AND PedigreeSId = ?",
                        undef, $MapId, $StudyId, $PedigreeSId);
            }
        }
    }
    (! defined($ped)) and error($KelvinDataset::errstr);
    
    if ($$study{role} eq "SERVER") {
        $dbh->do("CALL BadScaling(?)", undef, $StudyId);
        # Add the following to hold trait likelihood references
        $dbh->do("INSERT IGNORE INTO PedigreePositions "
                . "(StudyId, PedigreeSId, ChromosomeNo, RefTraitPosCM, "
                . "PedTraitPosCM, MarkerCount, FreeModels) "
                . "SELECT DISTINCT $StudyId, PedigreeSId, ChromosomeNo, "
                . "-9999.99, -9999.99, MarkerCount, 0 FROM PedigreePositions "
                . "WHERE StudyId = $StudyId", undef);
        return; # All client-only from here on out...
    }
    
    # 'Freshen' the Positions tables
    # Three cases we know of: marker, individual value, and range 
    # specification, and all can be in lists
    my $JointTPs = join(',', @{$config->isConfigured("TraitPositions")});
    $JointTPs =~ s/\s+//g;
    my @TPs = split(',',$JointTPs);
    $ChromosomeNo = $$dataset{chromosome};
    # First find where "begin" and "end" are for the current set of markers.
    my $lastMarkerPos = 0;
    my $firstMarkerPos = 9999.9;
    %markers = %{$$dataset{markers}};
    foreach my $marker (keys %markers) {
        $lastMarkerPos = $markers{$marker}{avgpos}
                if ($markers{$marker}{avgpos} > $lastMarkerPos);
        $firstMarkerPos = $markers{$marker}{avgpos}
                if ($markers{$marker}{avgpos} < $firstMarkerPos);
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
            $dbh->do("INSERT IGNORE INTO Positions (StudyId, ChromosomeNo, "
                    . "RefTraitPosCM) SELECT $StudyId, $ChromosomeNo, AvePosCM"
                    . " FROM MapMarkers WHERE StudyId = $StudyId "
                    . "AND MapId = $MapId AND AvePosCM NOT IN "
                    . "(SELECT DISTINCT RefTraitPosCM FROM Positions WHERE "
                    . "StudyId = $StudyId)", undef);
        } elsif ($TP =~ /(\d*.?\d*)-(\d*.?\d*):(\d*.?\d*)/) {
            my ($begin, $end, $inc, $precision, $va) = ($1, $2, $3, $3, 0);
            $precision =~ s/^[^\.]\.?(\d*)?/$1/;
            my $PosCM;
            do {
                $PosCM = sprintf("%0.*f", length($precision),
                        $begin + $inc * $va++);
                $dbh->do("INSERT IGNORE INTO Positions (StudyId, ChromosomeNo,"
                        . " RefTraitPosCM) VALUES (?,?,?)", undef, $StudyId,
                        $ChromosomeNo, $PosCM);
            } while ($PosCM < $end);
        } else {
            $dbh->do("INSERT IGNORE INTO Positions (StudyId, ChromosomeNo, "
                    . "RefTraitPosCM) VALUES (?,?,?)", undef, $StudyId,
                    $ChromosomeNo, $TP);
        }
    }
    
    # Random name for temp table
    my $temptable = "PPs_" . join ('',
            map { chr(int(rand() * 26) + 65) } (0 .. 7));
    
    # Turn off automatic error handling
    $dbh->{PrintError} = $dbh->{RaiseError} = 0;
    
    # Finally, freshen the PedigreePositions as needed...
    my $MPMarkers = ${$config->isConfigured("Multipoint")}[0];
    
    until ($dbh->do("CREATE TEMPORARY TABLE $temptable "
            . "SELECT a.StudyId, a.PedigreeSId, b.ChromosomeNo, "
            . "b.RefTraitPosCM, $MPMarkers 'MarkerCount' FROM Pedigrees a, "
            . "Positions b WHERE a.StudyId = $StudyId "
            . "AND a.StudyId = b.StudyId")) {
        check_mysql_retry($dbh);
    }
    until ($dbh->do("INSERT INTO PedigreePositions (StudyId, PedigreeSId, "
            . "ChromosomeNo, RefTraitPosCM, MarkerCount) "
            . "SELECT a.StudyId, a.PedigreeSId, a.ChromosomeNo, "
            . "a.RefTraitPosCM, a.MarkerCount FROM $temptable a "
            . "LEFT OUTER JOIN PedigreePositions b ON a.StudyId = $StudyId "
            . "AND a.StudyId = b.StudyId AND a.PedigreeSId = b.PedigreeSId AND"
            . " a.ChromosomeNo = b.ChromosomeNo AND "
            . "a.RefTraitPosCM = b.RefTraitPosCM WHERE b.StudyId IS NULL")) {
        check_mysql_retry($dbh);
    }
    until ($dbh->do("CALL BadScaling(?)", undef, $StudyId)) {
        check_mysql_retry($dbh);
    }
    
    return;
}

sub check_mysql_retry
{
    my ($handle) = @ARG;
    
    ($DBI::errstr =~ /try restarting transaction/)
            or die("DBI failure executing '", $handle->{Statement}, "'\n",
            "$DBI::errstr\n");
}

sub fatal
{
    die ("FATAL - ABORTING, @ARG");
}

sub error
{
    die ("ERROR - EXITING, @ARG\n");
}

sub warner
{
    warn ("WARNING, @ARG\n");
}

main() unless caller(0);
1;
