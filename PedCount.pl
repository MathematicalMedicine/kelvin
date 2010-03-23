#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(sum min);
use Data::Dumper;
use File::Basename;
use lib dirname($0);
use kf;

#####################################
#
# $Id$
#
# by Bill Valentine-Cooper, portions from work by John Burian
#
# Copyright 2008, Nationwide Children's Research Institute.  All
# rights reserved.  Permission is hereby given to use this software
# for non-profit educational purposes only.
#

$| = 1;    # Force flush of output as printed.

# Sanctioned globals
my $PairCount = 0;
my $LiabilityClasses = 1;
my $maf0           = 0; # True if we found missing minor allele frequencies
my $pedFile          = "pedpost.dat";
my $mapFile          = "mapfile.dat";
my $locusFile    = "datafile.dat";
my $frequencyFile      = "markers.dat";                 # Default input files

# Heavy-duty diagnostics
my $DiagSolos = 0;    # Writes out a separate file for each marker with actual copies of template pedigrees

# Command line option flags
my $help           = 0;
my $config         = 0;
my $pre            = 0;
my $post           = 0;
my $noparents      = 0;
my $XC             = 0;
my $imprinting     = 0;
my $liability      = 0;
my $bare           = 0;
my $count          = 0;
my $write          = "unspecified";
my $loops          = 0;
my $stats          = 0;
my $split          = 0;
my $subdirectories = 0;
my $nokelvin       = 0;
my $quiet          = 0;
my @include        = ();
my @exclude        = ();
my $WritePrefix    = "PC";

# Sorting routines
sub numerically { $a <=> $b }

sub numericIsh {
    if (($a . $b) =~ /^\d+$/) {
        $a <=> $b;
    } else {
        $a cmp $b;
    }
}

#####################################
#
sub dirPhenoCodes {
    $UnknownAffection = $Directives{PhenoCodes}[0];
    $Unaffected       = $Directives{PhenoCodes}[1];
    $Affected         = $Directives{PhenoCodes}[2];
}

#####################################
#
sub dirCountFile {
    die "Cannot generate counts for configuration that already has them."
      if ($count);
}

#####################################
#
sub dirImprinting {
    $imprinting = 1;
}

#####################################
#
sub dirLiabilityClasses {
    $liability        = 1;
    $LiabilityClasses = $Directives{LiabilityClasses}[0];
}

#####################################
#
sub dirQT {

    # If it's not already explicitly specified, set QT default for PhenoCodes
    if (!defined($Directives{PhenoCodes})) {
        $UnknownAffection = -99.99;
        $Unaffected       = -88.88;
        $Affected         = 88.88;
    }
}

#####################################
#
sub dirSexLinked {
    $XC = 1;
}

#####################################
#
# Evolved from a trios-only version by John Burian. See -help output for documentation.
#
sub bucketizePedigrees {

    # These are translations for bucket names generated from the algorithm
    # so if you change one, change both.

    my %NiceNames = (

        # It's alleles+affectation for CC
        '000000111' => 'ctrl11',
        '000000112' => 'case11',
        '000000121' => 'ctrl12',
        '000000122' => 'case12',
        '000000221' => 'ctrl22',
        '000000222' => 'case22',

        # It's alleles+affectation+sex XC CC
        '000100021111' => 'ctrl1',
        '000100021121' => 'case1',
        '000100021112' => 'ctrl11',
        '000100021122' => 'case11',
        '000100021212' => 'ctrl12',
        '000100021222' => 'case12',
        '000100022211' => 'ctrl2',
        '000100022221' => 'case2',
        '000100022212' => 'ctrl22',
        '000100022222' => 'case22',

        # It's John's encoding for trios
        '001001002' => 'T30',
        '001001112' => 'T18',
        '001001122' => 'T19',
        '001001222' => 'T20',
        '001111002' => 'T24',
        '001111112' => 'T05',
        '001111122' => 'T06',
        '001121002' => 'T27',
        '001121112' => 'T12',
        '001121122' => 'T13',
        '001121222' => 'T14',
        '001221002' => 'T29',
        '001221122' => 'T16',
        '001221222' => 'T17',
        '111111002' => 'T21',
        '111111112' => 'T01',
        '111121002' => 'T22',
        '111121112' => 'T02',
        '111121122' => 'T03',
        '111221002' => 'T23',
        '111221122' => 'T04',
        '121121002' => 'T25',
        '121121112' => 'T07',
        '121121122' => 'T08',
        '121121222' => 'T09',
        '121221002' => 'T26',
        '121221122' => 'T10',
        '121221222' => 'T11',
        '221221002' => 'T28',
        '221221222' => 'T15',
    );

    my $Type   = shift();    # Pedigree type for writing
    my $Prefix = shift();    # Uniqifying (what a word!) prefix for files

    # Verify that this is a 2pt analysis (default, so look for multipoint directives)
    die "Generation of counts not permitted for a multipoint analysis.\n"
      if (defined($Directives{MultiPoint}));

    # Verify that all markers are present and only biallelic...
    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        next if ($LociAttributes{$Name}{Type} eq "C");
        next if ($LociAttributes{$Name}{Type} eq "A");
        next if ($LociAttributes{$Name}{Type} eq "T");
        die "No allele information found for marker $Name for count generation.\n"
          if (!defined($LociAttributes{$Name}{Alleles}));
        die "Marker $Name not biallelic, not permitted for count generation.\n"
          if (scalar(@{ $LociAttributes{$Name}{Alleles}{OrderedList} }) > 2);
    }

    my %Buckets   = ();    # Fully-funkified bucket names with encoded everything
    my @Skippies  = ();    # Pedigrees copied on thru without bucketization (affects stats)
    my %Templates = ();    # Template pedigrees
    my $PedSeq    = 1;     # Template pedigree ID to keep them short

    # Look at each family...

    for my $Ped (sort numericIsh keys %Pedigrees) {

        my $memberCount = scalar(keys %{ $Pedigrees{$Ped} });

        # Qualify the family for inclusion in trio buckets by
        # verifying depth of 1
        my $SkipFlag = 0;
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (
                (
                    ($Dad ne $UnknownPerson) && (($Pedigrees{$Ped}{$Dad}{Dad} ne $UnknownPerson)
                        || ($Pedigrees{$Ped}{$Dad}{Mom} ne $UnknownPerson))
                )
                || ($Mom ne $UnknownPerson) && (($Pedigrees{$Ped}{$Mom}{Dad} ne $UnknownPerson)
                    || ($Pedigrees{$Ped}{$Mom}{Mom} ne $UnknownPerson))
              ) {
                push @Skippies, $Ped;
                print "Will copy multi-generation pedigree $Ped intact.\n";
                $SkipFlag = 1;
                last;
            }
        }
        next if ($SkipFlag);

        # Generate the family bucket for each marker pair
        for my $i (0 .. $PairCount - 1) {
            my @bucketList = ();

	    # A bucket name encoding is:
	    #
	    #    <Dad key> := <sorted Dad alleles>.<Dad aff>[.<Dad LC>][.<Dad Sex>]
	    #    <Mom key> := <sorted Mom alleles>.<Mom aff>[.<Mom LC>][.<Mom Sex>]
	    #    <Child key> := <sorted Child alleles>.<Child aff>[.<Child LC>][.<Child Sex>]
	    #    <bucket> := <lower parent key>.<higher parent key>.<Child 1 key>[.<Child 2 key>...]
	    #
	    my $BucketName = "";
            for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {

                my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
                # Skip the parents by checking for individual's Dad.
                if ($Dad ne $UnknownPerson) {
		    if ($BucketName eq "") {
			# Build the common parental portion from the first child's parent information
			my $DadKey = $Pedigrees{$Ped}{$Dad}{Mks}[$i];
			($DadKey eq '2 1') and $DadKey = '1 2';
			$DadKey .= $Pedigrees{$Ped}{$Dad}{Aff};
			$DadKey .= $Pedigrees{$Ped}{$Dad}{LC} if ($liability);
			$DadKey .= $Pedigrees{$Ped}{$Dad}{Sex} if
			    (defined($Directives{SexLinked}) || $XC || defined($Directives{Imprinting}) || $imprinting);
			my $Mom    = $Pedigrees{$Ped}{$Ind}{Mom};
			my $MomKey = $Pedigrees{$Ped}{$Mom}{Mks}[$i];
			($MomKey eq '2 1') and $MomKey = '1 2';
			$MomKey .= $Pedigrees{$Ped}{$Mom}{Aff};
			$MomKey .= $Pedigrees{$Ped}{$Mom}{LC} if ($liability);
			$MomKey .= $Pedigrees{$Ped}{$Mom}{Sex} if
			    (defined($Directives{SexLinked}) || $XC || defined($Directives{Imprinting}) || $imprinting);
			$BucketName = ($MomKey gt $DadKey) ? $DadKey . $MomKey : $MomKey . $DadKey;
		    }
		    # Build and add the non-parent individual (Child) key to what we have for the bucket name already
                    my $ChildKey = $Pedigrees{$Ped}{$Ind}{Mks}[$i];
                    ($ChildKey eq '2 1') and $ChildKey = '1 2';
                    $ChildKey .= $Pedigrees{$Ped}{$Ind}{Aff};
		    $ChildKey .= $Pedigrees{$Ped}{$Ind}{LC} if ($liability);
		    $ChildKey .= $Pedigrees{$Ped}{$Ind}{Sex} if
			(defined($Directives{SexLinked}) || $XC || defined($Directives{Imprinting}) || $imprinting);

                    $BucketName .= $ChildKey;
                }
            }
	    $BucketName =~ s/ //g;
            $Buckets{ $Loci[ $i ] . "_" . $BucketName }++;

#	    print "For marker ".$Loci [ $i ]." bucket ".$BucketName." gets pedigree ".sprintf("%003d\n", $Ped);
            if (!defined($Templates{$BucketName})) {
                $Templates{$BucketName}{Ped}    = $Ped;
                $Templates{$BucketName}{PedSeq} = sprintf("P%04d", $PedSeq++);
                $Templates{$BucketName}{PairID} = $i;
            }
        }
    }

    if (scalar(keys %Buckets) == 0) {
        print "Bucketizing cannot reduce your pedigree count.\n";
    } else {
        print sprintf(
            "Bucketizing reduces evaluation count from %d to %d, or by %2d%%\n",
            scalar(keys %Pedigrees) * $PairCount,
            (scalar(keys %Buckets) + (scalar(@Skippies) * $PairCount)),
            100 - (
                100 *
                  (scalar(keys %Buckets) + (scalar(@Skippies) * $PairCount)) /
                  (scalar(keys %Pedigrees) * $PairCount)
            )
        );
    }

    return if (!$write);

    # Now write-out at least the pedigree and counts

    if ($Type eq "POST") {
        open OUT, ">" . $Prefix . "Pedigrees.Dat";
	print OUT '# $Id$'; print OUT "\n";
    } else {
        open OUT, ">" . $Prefix . "Pedigrees.Pre";
    }

    # First the intact (skipped) pedigrees
    for my $Ped (@Skippies) {
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            print OUT sprintf("%4s %3s %3s %3s ", $Ped, $Ind, $Pedigrees{$Ped}{$Ind}{Dad}, $Pedigrees{$Ped}{$Ind}{Mom});
            print OUT sprintf("%3s %3s %3s ",
                $Pedigrees{$Ped}{$Ind}{Kid1},
                $Pedigrees{$Ped}{$Ind}{nPs},
                $Pedigrees{$Ped}{$Ind}{nMs})
              if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Aff} . "  ";
            if ($liability) {
                print OUT $Pedigrees{$Ped}{$Ind}{LC} . "  ";
            }
            print OUT join("  ", @{ $Pedigrees{$Ped}{$Ind}{Mks} })."   ";;
            print OUT "Ped: $Ped Per: $Ind\n";
        }
    }

    # Next the template pedigrees
    my $WarnAboutCaseControlFlag = 1; # Should we warn about a case/control run?
    for my $PB (sort numericIsh keys %Templates) {
        my $Ped      = $Templates{$PB}{Ped};
        my $PairID   = $Templates{$PB}{PairID};
        my $PedSeq   = $Templates{$PB}{PedSeq};
        my $NiceName = defined($NiceNames{$PB}) ? $NiceNames{$PB} : $PB;
        if (($NiceName =~ /case|ctrl/) && $WarnAboutCaseControlFlag) {
            $WarnAboutCaseControlFlag = 0;
            print "WARNING - Pedigree counting for case-control analysis is still under development!\n";
        }

        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            print OUT
              sprintf("%4s %3s %3s %3s ", $NiceName, $Ind, $Pedigrees{$Ped}{$Ind}{Dad}, $Pedigrees{$Ped}{$Ind}{Mom});
            print OUT sprintf("%3s %3s %3s ",
                $Pedigrees{$Ped}{$Ind}{Kid1},
                $Pedigrees{$Ped}{$Ind}{nPs},
                $Pedigrees{$Ped}{$Ind}{nMs})
              if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Aff} . "   ";
            if ($liability) {
                print OUT $Pedigrees{$Ped}{$Ind}{LC} . "   ";
            }
            my $Pair = $Pedigrees{$Ped}{$Ind}{Mks}[$PairID];
            my ($Left, $Right) = split(/ /, $Pair);
            for my $i (0 .. $PairCount - 1) {
                next if (!$LociAttributes{ $Loci[ $i ] }{Included});
                next if ($LociAttributes{ $Loci[ $i ] }{Type} =~ /^[AT]$/);
                if ($Left != 0) {
                    print OUT $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList}[$Left - 1] . " ";
                } else {
                    print OUT AttributeMissing . " ";
                }
                if ($Right != 0) {
                    print OUT $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList}[$Right - 1] . "  ";
                } else {
                    print OUT AttributeMissing . "  ";
                }
            }
            print OUT " Ped: $NiceName Per: $Ind\n";
        }
    }
    close OUT;

    # Diagnostic multiple template pedigrees! (Out-of-date)
    if ($DiagSolos) {
        for my $i (0 .. $PairCount - 1) {
            next if (!$LociAttributes{ $Loci[ $i ] }{Included});

            if ($Type eq "POST") {
                open OUT, ">" . $Prefix . "_solo_" . $Loci[ $i ] . "Pedigrees.Dat";
            } else {
                open OUT, ">" . $Prefix . "_solo_" . $Loci[ $i ] . "Pedigrees.Pre";
            }

            print OUT '# $Id$'; print OUT "\n";

            for my $PB (sort numericIsh keys %Templates) {
                my $FB = $Loci[ $i ] . "_" . $PB;
                if (!defined($Buckets{$FB})) {
                    next;
                } else {
                    my $Ped    = $Templates{$PB}{Ped};
                    my $PairID = $Templates{$PB}{PairID};
                    my $PedSeq = $Templates{$PB}{PedSeq};
                    print "Writing "
                      . $Buckets{$FB}
                      . " copies of "
                      . $PedSeq
                      . " for marker "
                      . $Loci[ $i ] . "\n";
                    for my $j (1 .. $Buckets{$FB}) {
                        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
                            print OUT sprintf(
                                "%4s %3s %3s %3s ",
                                $PB . "_" . $j,
                                $Ind,
                                $Pedigrees{$Ped}{$Ind}{Dad},
                                $Pedigrees{$Ped}{$Ind}{Mom}
                            );
                            print OUT sprintf("%3s %3s %3s ",
                                $Pedigrees{$Ped}{$Ind}{Kid1},
                                $Pedigrees{$Ped}{$Ind}{nPs},
                                $Pedigrees{$Ped}{$Ind}{nMs})
                              if ($Type eq "POST");
                            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
                            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
                            my $Pair = "  " . $Pedigrees{$Ped}{$Ind}{Mks}[$PairID];
                            print OUT $Pedigrees{$Ped}{$Ind}{Aff} . " "
                              . join(" ", $Pair x min($split, $PairCount)) . "   ";
                            print OUT "Ped: " . $PB . "_" . $j . " Per: $Ind\n";
                        }
                    }
                }
            }
            close OUT;
        }
    }

    if ($Type ne "POST") {
        system("makeped " . $Prefix . "Pedigrees.Pre " . $Prefix . "Pedigrees.Dat N");
    }

    # Finally the counts.
    open OUT, ">" . $Prefix . "Counts.Dat";

    print OUT '# $Id$'; print OUT "\n";

    print OUT "MARKER ";
    for my $PB (sort numericIsh keys %Templates) {
        my $NiceName = defined($NiceNames{$PB}) ? $NiceNames{$PB} : $PB;
        print OUT $NiceName . " ";
    }
    print OUT "\n";
    for my $i (0 .. $PairCount - 1) {
        next if (!$LociAttributes{ $Loci[ $i ] }{Included});
        print OUT $Loci[ $i ] . " ";
        for my $PB (sort numericIsh keys %Templates) {
            my $FB = $Loci[ $i ] . "_" . $PB;
            if (!defined($Buckets{$FB})) {
                print OUT "  0 ";
            } else {
                print OUT sprintf("%3d ", $Buckets{$FB});
            }
        }
        print OUT "\n";
    }
    close OUT;

    open OUT, ">" . $Prefix . "Data.Dat";
    print OUT "T Trait\n";
    print OUT "C liabilityClass\n" if ($liability);

    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
    }
    close OUT;

    # If there was no configuration, create all of the supporting files
    if ($config) {
        if ($maf0) {
            print "Encountered missing biallelic marker alleles, writing full marker data file.\n";

            open OUT, ">" . $Prefix . "Markers.Dat";

            print OUT '# $Id$'; print OUT "\n";

            for my $Name (@Loci) {
                next if (!$LociAttributes{$Name}{Included});
                if ($LociAttributes{$Name}{Type} eq "M") {
                    print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
                    print OUT "F";
		    for my $Allele (@{ $LociAttributes{$Name}{Alleles}{OrderedList} }) {
			print OUT sprintf(" %.8f", $LociAttributes{$Name}{Alleles}{$Allele}{Frequency});
                    }
                    print OUT "\n";
		}
	    }
        }
        close OUT;

        print "Remember to modify your configuration file to specify the new pedigree, companion and count files.\n"
          if (!$quiet);
        return;
    }

    open OUT, ">" . $Prefix . "Config.Dat";

    print OUT '# $Id$'; print OUT "\n";

    my $configPrefix = $Prefix;
    $configPrefix =~ s/.*\///;
    print OUT "PedigreeFile " . $configPrefix . "Pedigrees.Dat\n";
    print OUT "LocusFile " . $configPrefix . "Data.Dat\n";
    print OUT "FrequencyFile " . $configPrefix . "Markers.Dat\n";
    print OUT "MapFile " . $configPrefix . "Map.Dat\n";
    print OUT "CountFile " . $configPrefix . "Counts.Dat\n";
    print OUT "BayesRatioFile " . $configPrefix . "BR.Out\n";
    print OUT "PPLFile " . $configPrefix . "PPL.Out\n";

    if ($liability) {
	# Determine number of unique liability classes
	my @LiabilityList = ();
	for my $Ped (%Pedigrees) {
	    for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
		push @LiabilityList, $Pedigrees{$Ped}{$Ind}{LC};
	    }
	}
	my %Seen;
	my @Uniqed = grep !$Seen{$_}++, @LiabilityList;
	print OUT "LiabilityClasses ".scalar(@Uniqed)."\n" if ($liability);
    }


    print OUT <<EOF;

# The rest is the standard analysis grid...
FixedModels
Theta 0-0.5:0.05
DPrime -1-1:0.1
LD
DiseaseGeneFrequency 0.001, 0.999, 0.1-0.9:.1
Alpha 0.05-1.0:0.1
Penetrance DD 0.0-0.9:0.1, 0.999
Penetrance Dd 0.0-0.9:0.1, 0.999
Penetrance dd 0.0-0.9:0.1, 0.999
Constrain Penetrance DD >= Dd
Constrain Penetrance Dd >= dd
Constrain Penetrance DD != Dd, Dd != dd

EOF
    print OUT "SexLinked\n" if (defined($Directives{SexLinked}) || $XC);
    close OUT;

    open OUT, ">" . $Prefix . "Markers.Dat";

    print OUT '# $Id$'; print OUT "\n";

    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        if ($LociAttributes{$Name}{Type} eq "M") {
            print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
            print OUT "F";
	    for my $Allele (@{ $LociAttributes{$Name}{Alleles}{OrderedList} }) {
		print OUT sprintf(" %.8f", $LociAttributes{$Name}{Alleles}{$Allele}{Frequency});
	    }
            print OUT "\n";
        }	
    }
    close OUT;

    #
    open OUT, ">" . $Prefix . "Map.Dat";

    print OUT '# $Id$'; print OUT "\n";

    print OUT "CHR MARKER KOSAMBI\n";
    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        if ($LociAttributes{$Name}{Type} eq "M") {
            print OUT "1 $Name 1\n";
        }
    }
    close OUT;
}

#####################################
#
sub writeExpanded {

    my $Type   = shift();    # Pedigree type for writing
    my $Prefix = shift();    # Uniqifying (what a word!) prefix for files

    # Write-out the pedigree and counts

    if ($Type eq "POST") {
        open OUT, ">" . $Prefix . "Pedigrees.Dat";
    } else {
        open OUT, ">" . $Prefix . "Pedigrees.Pre";
    }

    print OUT '# $Id$'; print OUT "\n";

    for my $Ped (sort numericIsh keys %Pedigrees) {
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            print OUT sprintf("%4s %3s %3s %3s ", $Ped, $Ind, $Pedigrees{$Ped}{$Ind}{Dad}, $Pedigrees{$Ped}{$Ind}{Mom});
            print OUT sprintf("%3s %3s %3s ",
                $Pedigrees{$Ped}{$Ind}{Kid1},
                $Pedigrees{$Ped}{$Ind}{nPs},
                $Pedigrees{$Ped}{$Ind}{nMs})
              if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Aff} . "  ";
            if ($liability) {
                print OUT $Pedigrees{$Ped}{$Ind}{LC} . "  ";
            }
            my @Pairs = @{ $Pedigrees{$Ped}{$Ind}{Mks} }; # Pairs are of allele ordinals, not names
            for my $i (0 .. $PairCount - 1) {
                next if (!$LociAttributes{ $Loci[ $i ] }{Included});
                next if ($LociAttributes{ $Loci[ $i ] }{Type} =~ /^[AT]$/);
#                print OUT $Pairs[$i] . "  ";

                my ($Left, $Right) = split(/ /, $Pairs[$i]);

                print OUT $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList}[$Left - 1]." ";
                print OUT $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList}[$Right - 1]."  ";

#                my ($Left, $Right) = split(/ /, $Pairs[$i]);
#                if (defined($LociAttributes{ $Loci[ $i ] }{Alleles}{$Left})) {
#                    print OUT $LociAttributes{ $Loci[ $i ] }{Alleles}{$Left} . " ";
#                } else {
#                    print OUT AttributeMissing . " ";
#                }
#                if (defined($LociAttributes{ $Loci[ $i ] }{Alleles}{$Right})) {
#                    print OUT $LociAttributes{ $Loci[ $i ] }{Alleles}{$Right} . "  ";
#                } else {
#                    print OUT AttributeMissing . "  ";
#                }
            }
            print OUT "Ped: $Ped Per: $Ind\n";
        }
    }
    close OUT;

    if ($Type ne "POST") {
        system("makeped " . $Prefix . "_pedigrees.Pre " . $Prefix . "_pedigrees.Dat N");
    }

    open OUT, ">" . $Prefix . "Data.Dat";

    print OUT '# $Id$'; print OUT "\n";

    print OUT "T Trait\n";
    print OUT "C liabilityClass\n" if ($liability);
    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
    }
    close OUT;

    if ($config) {
        print "Remember to modify your configuration file to specify the new pedigree and companion files.\n";
        return;
    }

    open OUT, ">" . $Prefix . "Markers.Dat";

    print OUT '# $Id$'; print OUT "\n";

    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        if ($LociAttributes{$Name}{Type} eq "M") {
            print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
            print OUT "F";
	    for my $Allele (@{ $LociAttributes{$Name}{Alleles}{OrderedList} }) {
		print OUT sprintf(" %.8f", $LociAttributes{$Name}{Alleles}{$Allele}{Frequency});
            }
            print OUT "\n";
        }
    }
    close OUT;

    open OUT, ">" . $Prefix . "Map.Dat";

    print OUT '# $Id$'; print OUT "\n";

    print OUT "CHR MARKER KOSAMBI\n";
    for my $Name (@Loci) {
        next if (!$LociAttributes{$Name}{Included});
        if ($LociAttributes{$Name}{Type} eq "M") {
            print OUT "1 $Name 1\n";
        }
    }
    close OUT;
}

#####################################
# Process -include and -exclude options to indicate which markers are
# important.
#
sub doMarkerInclusion {

    # Not starting with a guarenteed clean slate, so make sure they're all on
    for my $Name (@Loci) {
        $LociAttributes{$Name}{Included} = 1 if ($LociAttributes{$Name}{Type} eq "M");
    }
    if (scalar(@include)) {

        # If we're doing inclusion, then turn them all off
        for my $Name (@Loci) {
            $LociAttributes{$Name}{Included} = 0 if ($LociAttributes{$Name}{Type} eq "M");
        }

        # Then turn on what is requested
        for my $Name (my @copy = @include) {
            $Name = sprintf("M%04d", $Name) if (!$config);
            if (defined($LociAttributes{$Name}{Included})) {
                $LociAttributes{$Name}{Included} = 1 if ($LociAttributes{$Name}{Type} eq "M");
            } else {
                print "Marker \"$Name\" specified in -include list not found!\n";
            }
        }
    }
    if (scalar(@exclude)) {

        # Turn off what is requested
        for my $Name (my @copy = @exclude) {
            $Name = sprintf("M%04d", $Name) if (!$config);
            if (defined($LociAttributes{$Name}{Included})) {
                $LociAttributes{$Name}{Included} = 0 if ($LociAttributes{$Name}{Type} eq "M");
            } else {
                print "Marker \"$Name\" specified in -exclude list not found!\n";
            }
        }
    }
}

#####################################
# Check for kelvin constraints. This is just a convenient front-end for
# kelvin as it should still check for absolute constraints, i.e. those
# that no-one should be violating. This first tier of checking is more
# to ensure that regular users don't try exotic analyses that might give
# misleading results or take millenia to complete.
#
sub kelvinLimits {

    # General limitations
    warn "Warning -- kelvin currently only supports biallelic disease models!\n"
      if (defined($Directives{DiseaseAlleles}) && ($Directives{DiseaseAlleles}[0] != 2));
}


#####################################
# Verify command line parameters, check flags and do what the user asks.
#
print '$Id$'; print "\n";

my $Help = <<EOF;

Leverages pedigree inheritance patterns to reduce computational complexity
of linkage analysis. Multiple pedigrees that conform to a single inheritance 
pattern are converted to a single representative pedigree and a weight
(count) for each marker. Currently works only for two-point analysis because
it would be real trick for there to be many pedigrees where multiple markers 
all fit the same inheritance pattern.

The traditional analyses for which this was designed are:
1. Case/Control, where unrelated genotyped and phenotyped individuals are
   divided into cases and controls and fall into one of 6 inheritance pattern
   buckets normally, 10 for X-chromosome analysis.
2. Trios are 3-member pedigrees with one affected child and both parents 
   treated as if unaffected. They fall into one of 30 inheritance pattern
    buckets normally, 44 for X-chromosome analysis.
3. ASPs (Affected Sibling Pairs) are essentially expanded trios, i.e. 4-member
   pedigrees with two affected children and both parents in any affectation
   category, which will fall first into one of the 3-member trio buckets for 
   the first child, and then a second 3-member bucket with the same parentage 
   for the second child.

All cases are handled by building a "bucket key" for each family by:

1. Building parent components by concatenating ordered marker genotype
   w/affectation status.
2. For each child concatenating ordered marker genotype w/affectation status
   and pushing onto list of children.
3. Ordering list of children and appending to ordered parental components.

X-chromosome analysis requires only that gender be included in the parental
and child components.

This handles all nuclear families. While enumerating every possible bucket for all
possible nuclear family combinations would be exhaustive, we use Perl hashes
to produce only the buckets needed. When we create a bucket, we keep track of
the last pedigree that fit into it so we can use that pedigree's genotypic and
phenotypic information to create the template pedigree for the bucket. We could
decode the bucket name and achieve the same goal, but this is easier and more
reliable.

EOF

my $Usage = <<EOF;

Usage "perl $0 [<flags>...] <input file>"

where <flags> are any of:

-help		Print a more philosophical explanation of purpose.
-config		The input file specified is a KELVIN configuration file. Otherwise
		it is assumed to be a pedigree file.
-pre		Pedigrees are in pre-MAKEPED format.
-post		Pedigrees are in post-MAKEPED format.
-noparents	Pedigrees are in pre-MAKEPED format with no columns for parents.
-XC		This is a sex-linked (X-chromosome) analysis (also "SexLinked" in configuration 
		file)
-imprinting	This in an imprinting analysis (also IMP in configuration file)
-liability	The pedigree file has a column after affection status for liability 
		class. If no KELVIN configuration file was specified, the template
                configuration file will include a count of unique liability classes
                seen on the LiabilityClasses directive.xs
-bare		The pedigree file has only individual, affection status and marker
		allele pairs columns.
-nokelvin	Skip verification that kelvin can handle the analysis.
-quiet		Suppress extraneous output (used by kelvin frontend)
-loops		Check for consanguinity and marriage loops and print them if found.
-stats		Print statistics on the make-up of the pedigree(s).
-count		Count genotypically identical pedigrees and print statistics.
-include=<list>	Process only the markers named in the list. For pedigree file-only
		runs, marker names are sequence numbers, e.g. 2,3,4,7 (no spaces) 
		and ranges like 2-4 can be specified as well. Can be specified
		multiple times and all will apply.
-exclude=<list> Process all markers except those named in the list. Can be specified
		multiple times and all will apply. If both -include and -exclude are
		specified, the processing order will be includes then excludes.
-split=<n>	Split a two-point multiple-marker analysis into subsets of
		<n> markers each. For multipoint analysis, subsets are by marker
		groups, i.e. all requested trait loci that use the same set of markers,
		and <n> is the number of marker groups per analysis.
-write=[<prefix>] Write new files, optionally with prefix <prefix> instead of "PC".
		If -count was specified, genotypically identical pedigrees will be
		reduced to a single representative and a count file will be produced. 
		If -split was specified, separate sequenced sets of files will be 
		produced.
-subdirectories	New files produced are written to subdirectories named with the prefix
		that would have otherwise been used for filenames, e.g. PC1/Pedigrees.Dat.

The input file will be read and analyzed. If it is a configuration file,
the input files specified in directives will be read and analyzed as well. 

Validates pedigrees and other configuration files, warns of affectation
inconsistencies, finds loops, adds parents to orphaned case-control 
individuals, generates count files for genotypically identical single-
generation pedigrees (case/control, trios, ASPs, etc.) Splits two-point
analyses into separate runs. Selects ranges of markers for inclusion or
exclusion.

If counts are requested with -write, a new pedigree file (in the same 
PRE/POST format as the input) will be generated along with a count file. 
If only a pedigree file was provided (no -config flag), then dummy configuration, 
pedigree companion, marker and map files will be generated for a two-point 
LD run with marker allele frequencies computed from the control individuals.
If a configuration file was provided, the marker names from the pedigree 
companion and marker files it specifies will be used in the count file, and
no other supporting files generated.

Output files produced by -write flag (<prefix> will be used instead of "PC" if 
specified):
PC1_Counts.Dat - counts of genotypically identical pedigrees.
PC1_Pedigrees.Pre or .Dat - pedigrees corresponding to counts, and any uncounted
        "pass-thru" pedigrees. If in pre-makeped format and makeped is found,
        it will be run to generate PC1_Pedigrees.Dat.

Output files generated by -write without the -config flag (<prefix> will be used
instead of "PC" if specified):
PC1_Data.Dat - pedigree companion data file indicating loci column usage.
PC1_Markers.Dat - marker data with allele freqs from unaffected (control) individuals.
PC1_Map.Dat - minimal dummy map file.
PC1_Config.Dat - template kelvin config file.

NOTE - If you intend to hand-combine the results of different counted runs, then 
you must provide a marker file in order to ensure that the same allele name 
mappings are made every time.

EOF

GetOptions(
    'help'           => \$help,
    'config'         => \$config,
    'pre'            => \$pre,
    'post'           => \$post,
    'noparents'      => \$noparents,
    'XC'             => \$XC,
    'imprinting'     => \$imprinting,
    'liability'      => \$liability,
    'bare'           => \$bare,
    'nokelvin'       => \$nokelvin,
    'quiet'          => \$quiet,
    'loops'          => \$loops,
    'stats'          => \$stats,
    'count'          => \$count,
    'include=s'      => \@include,
    'exclude=s'      => \@exclude,
    'split=i'        => \$split,
    'subdirectories' => \$subdirectories,
    'write:s'        => \$write,
) or die "Invalid command line parameters.";
if ($write ne "unspecified") {
    $WritePrefix = $write if ($write ne "");
    $write = 1;
} else {
    $write = 0;
}

$Data::Dumper::Sortkeys = 1;

print $Help if ($help);

die "Invalid number of arguments supplied.\n$Usage"       if (($#ARGV < 0) && (!$quiet));
print "-help flag seen\n"                                 if (($help) && (!$quiet));
print "-config flag seen\n"                               if (($config) && (!$quiet));
print "-pre flag seen\n"                                  if (($pre) && (!$quiet));
print "-post flag seen\n"                                 if (($post) && (!$quiet));
print "-noparents flag seen\n"                            if (($noparents) && (!$quiet));
print "-XC flag seen\n"                                   if (($XC) && (!$quiet));
print "-imprinting flag seen\n"                           if (($imprinting) && (!$quiet));
print "-liability flag seen\n"                            if (($liability) && (!$quiet));
print "-bare flag seen\n"                                 if (($bare) && (!$quiet));
print "-nokelvin flag seen\n"                             if (($nokelvin) && (!$quiet));
print "-loops flag seen\n"                                if (($loops) && (!$quiet));
print "-stats flag seen\n"                                if (($stats) && (!$quiet));
print "-count flag seen\n"                                if (($count) && (!$quiet));
print "-include list of " . Dumper(\@include) . " seen\n" if ((@include) && (!$quiet));
print "-exclude list of " . Dumper(\@exclude) . " seen\n" if ((@exclude) && (!$quiet));
print "-split of $split seen\n"                           if (($split) && (!$quiet));
print "-subdirectories flag seen\n"                       if (($subdirectories) && (!$quiet));
print "-write seen, using \"$WritePrefix\" prefix\n"      if (($write) && (!$quiet));
die "-pre -post and -bare are mutually exclusive flags."
  if ($pre + $post + $bare > 1);

# Setup the dispatch table for parsing configuration
$KnownDirectives{PhenoCodes} = \&dirPhenoCodes;
$KnownDirectives{CountFile} = \&dirCountFile;
$KnownDirectives{Imprinting}= \&dirImprinting;
$KnownDirectives{LiabilityClasses} = \&dirLiabilityClasses;
$KnownDirectives{QT} = \&dirQT;
$KnownDirectives{SexLinked} = \&dirSexLinked;

if ($config) {
    my $ConfFile = shift;
    loadConf($ConfFile);
    if (defined($Directives{LocusFile}[0])) { $locusFile = $Directives{LocusFile}[0]; }
    $liability = loadCompanion($locusFile);
    if (defined($Directives{MapFile}[0])) { $mapFile = $Directives{MapFile}[0]; }
    loadMap($mapFile);
    if (defined($Directives{FrequencyFile}[0])) { $frequencyFile = $Directives{FrequencyFile}[0]; }
    loadFrequencies($frequencyFile);
    if (defined($Directives{PedigreeFile}[0])) { $pedFile = $Directives{PedigreeFile}[0]; }
    $maf0 = addMissingAlleles();
} else {
    $pedFile = shift;
}
#print Dumper(\@Loci);
#print Dumper(\%LociAttributes);

my $Type = "";
if ($pre) { $Type = "PRE"; } elsif ($post) { $Type = "POST"; } elsif ($bare) { $Type = "BARE"; }
my $pedFileType = assessPedigree($pedFile, $Type, $liability, $noparents);
$PairCount = loadPedigree($pedFile, $pedFileType, $liability, $noparents);
#print Dumper(\%Pedigrees);

if (!$config) {
    # Make-up @Loci and %LociAttributes if we have to...
    deriveAlleleFrequencies();
    $maf0 = addMissingAlleles();
}

checkRelations($pedFileType);
checkIntegrity();

if (!$nokelvin) {
    kelvinLimits();
}

if ($stats) {
    perfStats();
}

if ($loops) {
    if (!consanguinityLoop()) {
        marriageLoop();
    }
}

@include = split(',', join(',', @include));    # No support for ranges
@exclude = split(',', join(',', @exclude));

#@include = expand(join(',', @include)) if (scalar(@include));    # Support for ranges
#@exclude = expand(join(',', @exclude)) if (scalar(@exclude));

doMarkerInclusion();

if (defined($Directives{Multipoint})) {
    if ($split) {
        die "The -split option is not yet implemented for multipoint.\n";
    } else {
        die "Generation of counts not permitted for a multipoint analysis.\n" if ($count);
        writeExpanded($pedFileType) if ($write);
    }
} else {

    # Set split to the count of marker loci if it's not specified
    # on the command line so we can "split" the analysis into one piece
    # and not have redundant code.
    $split = $PairCount if (!$split);
    my $IncludedMarkers = 0;    # Number of markers included thus-far
    my $SplitSet        = 0;    # Sequence number of set of markers
    for my $i (0 .. $PairCount - 1) {

        if (($LociAttributes{ $Loci[$i] }{Included}) && (++$IncludedMarkers >= $split)) {

            # Hit our limit, exclude all the rest
            for my $j (($i + 1) .. $PairCount - 1) {
                $LociAttributes{ $Loci[$j] }{Included} = 0;
            }
            # Do the work
            print "Marker set " . (++$SplitSet) . "\n";
            mkdir $WritePrefix . $SplitSet if ($write && $subdirectories);
            if ($count) {
                bucketizePedigrees($pedFileType, $WritePrefix . $SplitSet . (($write && $subdirectories) ? "/" : "_"));
            } elsif ($write) {
                writeExpanded($pedFileType, $WritePrefix . $SplitSet . (($write && $subdirectories) ? "/" : "_"));
            }

            # Redo the inclusion...
            doMarkerInclusion();

            # ...and exclusion up to current
            for my $j (0 .. $i) {
                $LociAttributes{ $Loci[$j] }{Included} = 0;
            }
            $IncludedMarkers = 0;
        }
    }
    if ($IncludedMarkers != 0) {

        # Do the rest
        mkdir $WritePrefix . $SplitSet if ($write && $subdirectories);
        print "Marker set " . (++$SplitSet) . "\n";
        if ($count) {
            bucketizePedigrees($pedFileType, $WritePrefix . $SplitSet . (($write && $subdirectories) ? "/" : "_"));
        } elsif ($write) {
            writeExpanded($pedFileType, $WritePrefix . $SplitSet . (($write && $subdirectories) ? "/" : "_"));
        }
    }
}
exit;
