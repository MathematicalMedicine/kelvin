#!/usr/bin/env perl

# Creates a new Study within a LKS database and adds that which the study will
# be analyzing to same.
# 
# Author: Jo Valentine-Cooper <jvc@mathmed.org>
# Based on InitStudy.pl by Bill Valentine-Cooper
# 
# Copyright (C) 2015, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

package LKS_InitStudy v2.0.1;

use v5.20;
use warnings 'FATAL';
use lib ($ENV{TOOLPATH} ? split (/:/, $ENV{TOOLPATH}) : ' ');
use English;

use DBI; # Database interaction
use DBI qw(:sql_types);

use File::Basename; # Get map file name and not the 128+ character path
use Digest::MD5;
use IO::File;

use BCMM::CLIParser;
use BCMMTools;
use RunSQLScript v1.1.2;
use KelvinDataset v1.4.0;
use KelvinConfig v1.5.0;
use KelvinFamily v1.4.0;

use Exporter qw(import);
our @EXPORT_OK = qw(init_study_tables insert_initial_study_data);

$OUTPUT_AUTOFLUSH= 1 ;  # Show the output when I say so.
my $svnid = '$Id$';

my $optssetup = [
    "Creates a new Study within a LKS database and adds that which the study "
            . "will be analyzing to same.",
    
    {   # custom option --clientconf (client config file)
        NAME => "clientconf",
        REQUIRED => 1,
        ARGTYPE => "client config file",
        HELP => "Specifies the Kelvin LKS client config file to use",
        VALIDATE => \&BCMM::CLIParser::validate_file_arg
    },
    {   # custom option --serverconfs (server config files)
        NAME => "serverconfs",
        REQUIRED => 1,
        ARGTYPE => "serverconf,serverconf,...",
        HELP => "Comma-separated list of Kelvin LKS server configs to use",
        VALIDATE => (
            sub {
                BCMM::CLIParser::validate_file_arg(@ARG, undef, "listmode");
            }
        ),
    },
    "outprefix",
    {   # custom option --checkoutfile (checksum output file)
        NAME => "checkoutfile",
        ARGTYPE => "checksums output file",
        HELP => "File to output checksums of our conf files to, to serve as a "
                . "touchstone for prereq-checking systems like depcheck",
    }
];


sub main {
    $main::VERSION = our $VERSION;  # needed for --version option to work
    init_study(\@ARGV);
}

sub init_study {
    my $opts = BCMM::CLIParser::autoparse_cmdline($optssetup, \@ARG);
    
    # potentially prep for checksum output
    my $checkout;
    if ($$opts{checkoutfile}) {
        $checkout = IO::File->new(
                "$$opts{outprefix}/$$opts{checkoutfile}.working", "w")
                or error("could not open checksum output file: $!");
    }
    
    foreach my $conffile (($$opts{clientconf}, @{$$opts{serverconfs}})) {
        my $config = KelvinConfig->new($conffile)
                or error("new KelvinConfig failed: $KelvinConfig::errstr");
        
        $config->validate or error("$KelvinConfig::errstr");
        $config->isConfigured("Study")
                or error("Study directive must be provided");
        $config->isConfigured("Multipoint")
                or error("Must be a multipoint analysis");
        
        my $study = $config->readStudyLine();
        if ($$study{role} eq "CLIENT") {
            init_study_tables($$study{host}, $$study{database},
                    $$study{username}, $$study{password}, $$opts{verbose});
        }
        insert_initial_study_data($config);
        
        if (defined($checkout)) {
            my $checksum = Digest::MD5->new();
            my $confdata = IO::File->new($conffile, "r")
                    or error("could not open $conffile: $!");
            $checksum->addfile($confdata);
            $checkout->print("$conffile: " . $checksum->hexdigest . "\n");
            $confdata->close();
        }
    }
    
    if (defined($checkout)) { $checkout->close(); }
    rename("$$opts{outprefix}/$$opts{checkoutfile}.working",
            "$$opts{outprefix}/$$opts{checkoutfile}");
}

sub init_study_tables {
    # Given the login info for a database, sets up that database with the
    # tables and procedures for LKS.
    
    my ($host, $dbname, $user, $pass, $verbosemode) = @ARG;
    
    # SQL scripts that need to be run, in this order, to initialize a new LKS
    # database
    my $sqlscriptnames = [
        ["LKS_setup_tables.sql", undef],
        ["LKS_setup_trigger_proc.sql", undef],
        ["DModelParts.sql", "lazyparse"],
        ["QModelParts.sql", "lazyparse"],
    ];
    
    my $dbh = RunSQLScript::connect_to_db($host, $dbname, $user, $pass,
            undef, undef, undef, $verbosemode);
    my $rowcount = $dbh->do("Set SQL_MODE=STRICT_ALL_TABLES");
    
    foreach my $sqlscriptpair (@$sqlscriptnames) {
        my ($sqlscriptfile, $lazyparse) = @$sqlscriptpair;
        my $sqlscript = BCMMTools::find_in_path($sqlscriptfile);
        my $sqlcmds = RunSQLScript::split_sql_script($sqlscript, $lazyparse);
        RunSQLScript::run_sql_cmdlist($dbh, $sqlcmds, undef, $verbosemode);
    }
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

sub insert_initial_study_data {
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
    
    # If we're running a client, clean out our study just for safety's sake
    # (we wait until now just in case something goes wrong with setting up the
    # KelvinDataset - don't want to touch the DB until we're sure we have what
    # we need)
    if ($$study{role} eq "CLIENT") {
        $dbh->do("CALL CleanStudy('$$study{label}')");
    }
    
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
    $dbh->do("CALL GetMapId(\'$StudyId\', \'"
            . (($$study{role} eq "CLIENT") ? 1 : 0)
            . "\', \'$MapScale\', \'".basename($MapFile)."\', \@MapId, \@MapScale, "
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
    } else {
        # 'Freshen' the Positions tables
        my $TPs = $dataset->expand_traitpositions(
                ${$config->isConfigured("TraitPositions")}[0]);
        $ChromosomeNo = $$dataset{chromosome};
        foreach my $TP (@$TPs) {
            $dbh->do("INSERT IGNORE INTO Positions (StudyId, ChromosomeNo, "
                    . "RefTraitPosCM) VALUES (?,?,?)", undef, $StudyId,
                    $ChromosomeNo, $TP);
        }
        
        # Random name for temp table
        my $temptable = "PPs_" . join ('',
                map { chr(int(rand() * 26) + 65) } (0 .. 7));
        
        # Turn off automatic error handling for these final statements
        $dbh->{PrintError} = $dbh->{RaiseError} = 0;
        
        # Finally, freshen the PedigreePositions as needed...
        my $MPMarkers = ${$config->isConfigured("Multipoint")}[0];
        
        until ($dbh->do("CREATE TEMPORARY TABLE $temptable "
                . "SELECT a.StudyId, a.PedigreeSId, b.ChromosomeNo, "
                . "b.RefTraitPosCM, $MPMarkers 'MarkerCount' FROM "
                . "Pedigrees a, Positions b WHERE a.StudyId = $StudyId "
                . "AND a.StudyId = b.StudyId")) {
            check_mysql_retry($dbh);
        }
        until ($dbh->do("INSERT INTO PedigreePositions (StudyId, PedigreeSId, "
                . "ChromosomeNo, RefTraitPosCM, MarkerCount) "
                . "SELECT a.StudyId, a.PedigreeSId, a.ChromosomeNo, "
                . "a.RefTraitPosCM, a.MarkerCount FROM $temptable a "
                . "LEFT OUTER JOIN PedigreePositions b "
                . "ON a.StudyId = $StudyId AND a.StudyId = b.StudyId "
                . "AND a.PedigreeSId = b.PedigreeSId "
                . "AND a.ChromosomeNo = b.ChromosomeNo "
                . "AND a.RefTraitPosCM = b.RefTraitPosCM "
                . "WHERE b.StudyId IS NULL")) {
            check_mysql_retry($dbh);
        }
        until ($dbh->do("CALL BadScaling(?)", undef, $StudyId)) {
            check_mysql_retry($dbh);
        }
        
    }
}

sub check_mysql_retry
{
    my ($handle) = @ARG;
    
    ($DBI::errstr =~ /try restarting transaction/)
            or die("DBI failure executing '", $handle->{Statement}, "'\n",
            "$DBI::errstr\n");
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
