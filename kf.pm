#!/usr/bin/perl -w
#use strict;

package kf;

use List::Util qw(sum min);
use Data::Dumper;
require Exporter;
our @ISA = qw(Exporter);

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
# @todo - Deal with actual column ordering in so-called "locus" file.
#

$| = 1;    # Force flush of output as printed.

# Permanent defaults
use constant AttributeMissing => "0";    # For marker alleles and Sex

# Sanctioned globals

my $HaveConfig         = 0;
my $XC             = 0;

# Defaults to be overridden by configuration file directives
our $UnknownAffection = 0;
our $Unaffected       = 1;
our $Affected         = 2;                             # Default affection indicators
our $UnknownPerson    = "0";                           # Default unknown person indicator
our $MapFunction      = "Kosambi, bless his heart!";

# Read/calculated results
our %Pedigrees;                                        # Pedigrees as loaded
our %Directives;                                       # Directives as loaded
our $PairCount      = 0;                                                # Last pedigree count of marker pairs
our @Loci;                                        # Ordered loci name list from companion file
our %LociAttributes;  # Loci attributes from companion and marker files
our %Map;                                                               # Loci on the map and other map attributes

# Nuisances to fix
my @Depths;                                                            # Referenced recursively and I'm fuddled
my @Ancestors;                                                         # ditto
my $ShortestLoop = "";    # Just batted around too much to monkey with right now.

# Known directives as well as dispatch routine if needed. I could avoid the NoActions, but
# I prefer to be explicit.
our %KnownDirectives = (
			"FrequencyFile" => \&NoAction,
			"MapFile" => \&NoAction,
			"PedigreeFile" => \&NoAction,
			"LocusFile" => \&NoAction,
			"BayesRatioFile" => \&NoAction,
			"PPLFile" => \&NoAction,
			"CountFile" => \&NoAction,
			"MODFile" => \&NoAction,
			"NIDetailFile" => \&NoAction,

			"NonPolynomial" => \&NoAction,
			"Imprinting" => \&NoAction,
			"SexLinked" => \&NoAction,
			"FixedModels" => \&NoAction,
			"DryRun" => \&NoAction,
			"ExtraMODs" => \&NoAction,

			"PolynomialScale" => \&NoAction,
			"LiabilityClasses" => \&NoAction,
			"DiseaseAlleles" => \&NoAction,

			"TraitPositions" => \&NoAction,
			"DiseaseGeneFrequency" => \&NoAction,
			"DPrime" => \&NoAction,
			"Theta" => \&NoAction,
			"MaleTheta" => \&NoAction,
			"FemaleTheta" => \&NoAction,
			"Alpha" => \&NoAction,
			"Penetrance" => \&NoAction,
			"Constrain" => \&NoAction,
			"Multipoint" => \&NoAction,
			"MarkerToMarker" => \&NoAction,
			"SexSpecific" => \&NoAction,
			"LD" => \&NoAction,
			"QT" => \&NoAction,
			"QTT" => \&NoAction,
			"Mean" => \&NoAction,
			"StandardDev" => \&NoAction,
			"DegOfFreedom" => \&NoAction,
			"Threshold" => \&NoAction,
			"Truncate" => \&NoAction,
			"PhenoCodes" => \&NoAction,
			"SurfacesPath" => \&NoAction,
			"Log" => \&NoAction,
			);

our @EXPORT = qw(
		 %Map @Loci %LociAttributes %Pedigrees %Directives %KnownDirectives 
		 AttributeMissing $UnknownAffection $UnknownPerson $Affected $Unaffected
		 addMissingAlleles assessPedigree deriveAlleleFrequencies checkRelations 
		 checkIntegrity loadPedigree loadConf loadCompanion loadFrequencies loadMap
		 );

#####################################
#
sub NoAction {
}

#####################################
#
sub followPaths {
    my ($Pedily, $Start, $End, $Seen, $i, $Ind);
    $Pedily = shift();
    $Start  = shift();
    $End    = shift();
    $Seen   = shift();

    # Have we ever reached our goal?
    if ($ShortestLoop ne "") {
        if (length($ShortestLoop) <= length($Seen)) {

#	    print "Overly-long path\n";
            return;
        }
    }

    # Have we reached our goal this time?
    if ($Start eq $End) {

        # Woohoo!
        $Seen .= "+" . $Start;
        if (!$ShortestLoop ne "") {
            $ShortestLoop = $Seen;
        } else {
            if (length($ShortestLoop) > length($Seen)) {
                $ShortestLoop = $Seen;
            }
        }
        return;
    }

    # Have we seen this one before?
    if ($Seen =~ /$Start/) {

        # Dead end
#	print "Dead end [$Seen]\n";
        return;
    }

    $Seen .= "+" . $Start;

    # Try up
    my ($Dad, $Mom);
    $Dad = $Pedigrees{$Pedily}{$Start}{Dad};
    if ($Dad != $UnknownPerson) {
        followPaths($Pedily, $Dad, $End, $Seen);
    }
    $Mom = $Pedigrees{$Pedily}{$Start}{Mom};
    if ($Mom != $UnknownPerson) {
        followPaths($Pedily, $Mom, $End, $Seen);
    }

    # Try down (for every individual who has $Start as a parent)
    for $Ind (keys %{ $Pedigrees{$Pedily} }) {
        if (   ($Pedigrees{$Pedily}{$Ind}{Dad} eq $Start)
            || ($Pedigrees{$Pedily}{$Ind}{Mom} eq $Start)) {
            followPaths($Pedily, $Ind, $End, $Seen);
        }
    }

    # Try my spouses (for every individual I share parenting with)
    my @Spouses    = ();
    my %SpouseSeen = ();
    for $Ind (keys %{ $Pedigrees{$Pedily} }) {
        if ($Start eq $Pedigrees{$Pedily}{$Ind}{Mom}) {
            push @Spouses, $Pedigrees{$Pedily}{$Ind}{Dad} unless $SpouseSeen{ $Pedigrees{$Pedily}{$Ind}{Dad} }++;
        }
        if ($Start eq $Pedigrees{$Pedily}{$Ind}{Dad}) {
            push @Spouses, $Pedigrees{$Pedily}{$Ind}{Mom} unless $SpouseSeen{ $Pedigrees{$Pedily}{$Ind}{Mom} }++;
        }
    }
    for my $Spouse (@Spouses) {
        followPaths($Pedily, $Spouse, $End, $Seen);
    }

    # Try over (for every individual with a same parent as me but me)
    # but if both parents are the same and the sibling has no kids, don't bother.
    for $Ind (keys %{ $Pedigrees{$Pedily} }) {
        my $YourDad = $Pedigrees{$Pedily}{$Ind}{Dad};
        my $MyDad   = $Pedigrees{$Pedily}{$Start}{Dad};
        my $YourMom = $Pedigrees{$Pedily}{$Ind}{Mom};
        my $MyMom   = $Pedigrees{$Pedily}{$Start}{Mom};
        if (   (($YourMom ne $UnknownPerson) && ($YourMom eq $MyMom))
            || (($YourDad ne $UnknownPerson) && ($YourDad eq $MyDad))) {
            if (($YourMom eq $MyMom) && ($YourDad eq $MyDad)) {
                for my $Child (keys %{ $Pedigrees{$Pedily} }) {
                    if (   ($Pedigrees{$Pedily}{$Child}{Dad} eq $Ind)
                        || ($Pedigrees{$Pedily}{$Child}{Mom} eq $Ind)) {
                        followPaths($Pedily, $Ind, $End, $Seen);
                        last;
                    }
                }
            } else {
                followPaths($Pedily, $Ind, $End, $Seen);
            }
        }
    }

#    print "End of road [$Seen]\n";
    return;
}

#####################################
# Recursively populate @Depths and @Ancestors
#
sub listAncestors {
    my ($Ped, $Ind, $Depth);
    $Ped   = shift();
    $Ind   = shift();
    $Depth = shift();

    $Depth++;
    my ($Dad, $Mom);
    $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
    $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
    if ($Dad ne $UnknownPerson) {
        push @Depths,    $Depth;
        push @Ancestors, $Dad;
        listAncestors($Ped, $Dad, $Depth);
    }
    if ($Mom ne $UnknownPerson) {
        push @Depths,    $Depth;
        push @Ancestors, $Mom;
        listAncestors($Ped, $Mom, $Depth);
    }
    return;
}

#####################################
#
sub countAncestors {
    my ($Pedily, $Individual, $Count);
    $Pedily     = shift();
    $Individual = shift();
    $Count      = shift();

    my ($Dad, $Mom, $DadCount, $MomCount);
    $Dad = $Pedigrees{$Pedily}{$Individual}{Dad};
    $Mom = $Pedigrees{$Pedily}{$Individual}{Mom};
    if ($Dad eq $UnknownPerson) {
        $DadCount = 0;
    } else {
        $DadCount = 1 + countAncestors($Pedily, $Dad, $Count);
    }
    if ($Mom eq $UnknownPerson) {
        $MomCount = 0;
    } else {
        $MomCount = 1 + countAncestors($Pedily, $Mom, $Count);
    }
    return ($DadCount + $MomCount + $Count);
}

#####################################
# One type of loop is two parents of a child, where the parents share common
# ancestry. Loop span is sum of two distances to common ancestor, so brother
# and sister parents would have a span of 2, because they both count back 1 to
# their mother (or father).
#
# Essentially we get all individuals in the ancestry of each pair and if
# there's commonality, then there's a loop.
#
sub consanguinityLoop() {
    my $LoopCount = 0;
    for my $Ped (keys %Pedigrees) {
        my @Pairings = ();
        my %Seen     = ();
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom ne $UnknownPerson) && ($Dad ne $UnknownPerson)) {
                push @Pairings, $Mom . "+" . $Dad unless $Seen{ $Mom . "+" . $Dad }++;
            }
        }
        for my $Pair (@Pairings) {
            my ($Mom, $Dad) = split /\+/, $Pair;
            @Depths    = ();
            @Ancestors = ();    # This could be better, but I'm cowardly
            listAncestors($Ped, $Mom, 0);
            my @MomDepths    = @Depths;
            my @MomAncestors = @Ancestors;
            @Depths    = ();
            @Ancestors = ();
            listAncestors($Ped, $Dad, 0);
            my @DadDepths    = @Depths;
            my @DadAncestors = @Ancestors;

            for my $i (0 .. $#MomDepths) {
                for my $j (0 .. $#DadDepths) {
                    if ($MomAncestors[$i] eq $DadAncestors[$j]) {
                        my $LoopSize = $MomDepths[$i] + $DadDepths[$j];
                        print "Pedigree $Ped consanguinity loop of size $LoopSize at ancestor "
                          . $MomAncestors[$i]
                          . " for pair $Mom / $Dad\n";
                        $LoopCount++;
                    }
                }
            }
        }
    }
    return $LoopCount;
}

#####################################
sub marriageLoop() {

# The more extensive loop that cares about more than common ancestry is any path
# over or up from one parent back down or over to the other without retracing steps.
# This approach can get a little slow at times.
#
    for my $Ped (keys %Pedigrees) {
        my %Seen     = ();
        my @Pairings = ();
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom ne $UnknownPerson) && ($Dad ne $UnknownPerson)) {
                push @Pairings, $Mom . "+" . $Dad unless $Seen{ $Mom . "+" . $Dad }++;
            }
        }
        for my $Pair (@Pairings) {
            my ($Mom, $Dad) = split /\+/, $Pair;

            # Follow paths from Mom's Dad looking for Dad, having seen Mom
            my $DadsGrandPa = $Pedigrees{$Ped}{$Dad}{Dad};
            my $DadsGrandMa = $Pedigrees{$Ped}{$Dad}{Mom};
            if (   ($DadsGrandPa ne $UnknownPerson)
                || ($DadsGrandMa ne $UnknownPerson)) {

                my $MomsGrandPa = $Pedigrees{$Ped}{$Mom}{Dad};
                my $MomsGrandMa = $Pedigrees{$Ped}{$Mom}{Mom};
                if ($MomsGrandPa ne $UnknownPerson) {
                    followPaths($Ped, $MomsGrandPa, $Dad, $Mom);
                }

                # Follow paths from Mom's Mom looking for Dad, having seen Mom
                if ($MomsGrandMa ne $UnknownPerson) {
                    followPaths($Ped, $MomsGrandMa, $Dad, $Mom);
                }
                if ($ShortestLoop ne "") {
                    print "Pedigree $Ped marriage loop for pair $Pair as [$ShortestLoop]\n";
                }
                $ShortestLoop = "";
            }
        }
    }
}

#####################################
# Open and read a kelvin configuration file producing a hash by directive of
# parameters in an array.
#
sub loadConf {
    $HaveConfig = 1;
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open $File\n";
    my $LineNo = 0;
    %Directives = ();
    while (<IN>) {
	$LineNo++;
	s/\s*\#.*//g;    # Trim comments before looking for semis.
	for (split /;/) {    # Semi is a directive delimiter -- not as good as a newline
	    s/^\s*//g;       # Trim leading whitespace
	    next if (/^$/);  # Drop empty lines
	    my @Parameters = split /[,\s]+/;
	    my $Directive  = shift @Parameters;
	    die "Configuration line $LineNo: \"$Directive\" is not a known configuration directive."
		if (!defined($KnownDirectives{$Directive}));
	    $Directives{$Directive} = \@Parameters;
	    &{ $KnownDirectives{$Directive} };
	}
    }
}

#####################################
# Open and assess a pedigree file to determine processing required. Return a type, one of
# POST, PRE or BARE.
#
sub assessPedigree {
    my $File = shift();
    my $Type = shift();
    my $liability = shift();
    my $noparents = shift();
	
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    while (<IN>) {
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        close IN;
        my @Columns      = split /\s+/;
        my $TotalColumns = scalar(@Columns);

        # See if the number of columns makes any kind of sense...
        my $MarkerColumns = (scalar(@Loci) - 1) * 2;
        if ($Type eq "PRE") {
            $TotalColumns += 2 if ($noparents);    # Make up for no parents here
            $TotalColumns -= 1 if ($liability);    # Make up for extra liability class column
                  # pre-MAKEPED is Ped Ind Dad Mom Sex Aff Pairs, i.e. some even number > 6
            die "Invalid number of columns in pre-MAKEPED pedigree file $File."
              if (($TotalColumns < 8) || ($TotalColumns & 1));
            die
              "Inconsistent column counts between companion datafile ($MarkerColumns markers) and pre-MAKEPED pedigree file $File ("
              . ($TotalColumns - 6) . ")."
              if ($HaveConfig && (($TotalColumns - $MarkerColumns) != 6));
            return "PRE";
        } elsif ($Type eq "POST") {
            $TotalColumns += 2 if ($noparents);    # Make up for no parents here
            $TotalColumns -= 1 if ($liability);    # Make up for extra liability class column
                  # post-MAKEPED is Ped Ind Dad Mom Kid1 nPs nMs Sex Prb Aff Pairs, i.e. some even number > 10
            die "Invalid number of columns ($TotalColumns) in post-MAKEPED pedigree file $File."
              if (($TotalColumns < 12)); # || ($TotalColumns & 1));
            die
              "Inconsistent column counts ($MarkerColumns of $TotalColumns total columns should be markers) between companion datafile and post-MAKEPED pedigree file $File."
              if ($HaveConfig && (($TotalColumns - $MarkerColumns) != 10) && (($TotalColumns - $MarkerColumns) != 14));
            return "POST";
        } elsif ($Type eq "BARE") {

            # Bare individuals for CC are $Ped/$Ind Aff Pairs, i.e. some even number > 2
            die "Invalid number of columns in bare pedigree file $File."
              if (($TotalColumns < 4) || ($TotalColumns & 1));
            die
              "Inconsistent column counts between companion datafile and post-MAKEPED pedigree file $File."
              if ($HaveConfig && (($TotalColumns - $MarkerColumns) != 2));
            return "BARE";
        } elsif ($MarkerColumns != 0) {
            $TotalColumns += 2 if ($noparents);    # Make up for no parents here
            my $Slack = $TotalColumns - $MarkerColumns;
            return "POST" if (($Slack == 14) || ($Slack == 10));    # With or without old "Ped: x Ind: y"
            return "PRE"  if ($Slack == 1);
            return "BARE" if ($Slack == 1);

            # Assume post-MAKEPED and be surprised if it isn't.
            return "POST";
        } else {

            # Assume post-MAKEPED and be surprised if it isn't.
            print
              "Unable to determine pedigree type based solely upon column count $TotalColumns, assuming post-MAKEPED.\n";
            return "POST";
        }
    }
}

#####################################
# Open and read a post-makeped pedigree file producing the three-dimensional hash
# %Pedigrees with the first being the family (pedigree) ID, the second being the
# individual ID, and the third being attributes of that individual.
#
sub loadPedigree {
    my $File = shift();
    my $Type = shift();
    my $liability = shift();
    my $noparents = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo      = 0;
    my $ColumnCount = 0;
    %Pedigrees = ();
    my $AlC = 0;    # Number of alleles
    my $GtC = 0;    # Number of genotyped alleles
    my ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, $LC, @Alleles);

    while (<IN>) {
        $LineNo++;
        print "Read $LineNo lines of $File\n" if (($LineNo % 1001) == 1000);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        my @Columns = split /\s+/;
        $ColumnCount = scalar(@Columns) if (!$ColumnCount);
        die "Inconsistent number of columns in $File at line $LineNo, was $ColumnCount, is " . scalar(@Columns) . "\n"
          if ($ColumnCount != scalar(@Columns));

        if ($Type eq "POST") {
            if ($noparents) {
                ($Ped, $Ind, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Alleles) = @Columns;
            } else {
                ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Alleles) = @Columns;
            }
            die "Pedigree line $LineNo: proband must be numeric, not \"$Prb\"." if (!($Prb =~ /[0-9]/));
        } elsif ($Type eq "PRE") {
            if ($noparents) {
                ($Ped, $Ind, $Sex, $Aff, @Alleles) = @Columns;
            } else {
                ($Ped, $Ind, $Dad, $Mom, $Sex, $Aff, @Alleles) = @Columns;
            }
        } elsif ($Type eq "BARE") {
            ($Ind, $Aff, @Alleles) = @Columns;
            $Ped = $Ind;
        } else {
            die "Unhandled pedigree type \"$Type\".";
        }

        # Validate everything we've got so far
        die "Pedigree line $LineNo: sex must be 1 or 2, not \"$Sex\"." if (($Type ne "BARE") && !($Sex =~ /[12]/));
        die
          "Pedigree line $LineNo: affection status must be $UnknownAffection, $Unaffected, or $Affected, not \"$Aff\"."
          if ( (!defined($Directives{QT}))
            && (($Aff != $UnknownAffection) && ($Aff != $Unaffected) && ($Aff != $Affected)));

        if ($liability) {
            $LC = shift @Alleles;
	    $Pedigrees{$Ped}{$Ind}{LC} = $LC;
        }

        # Handle marker pairs
        my ($OldFam, $OldInd, $Left);
        $AlC = 0;
        $GtC = 0;
        my @Pairs = ();
        for my $Allele (my @copy = @Alleles) {
            if ($Allele eq "Ped:") {
                $OldFam = $Alleles[ $AlC + 1 ];
                $OldInd = $Alleles[ $AlC + 3 ];
                last;
            }
            $AlC++;
	    my $Offset = (($AlC & 1) ? $AlC+1 : $AlC) / 2 - 1;
            my $Name;
            if (!$HaveConfig) {
                $Name = sprintf("M%04d", $Offset + 1);    # Offset == Name
                push @Loci, "$Name" if (scalar(@Loci) <= $Offset);
            } else {
                $Name = $Loci[ $Offset ];           # Integer division, thank you
            }
            if ($Allele ne AttributeMissing) {
                $GtC++;    # Keep track of how much genotypic information we have for this individual
                if (!$HaveConfig) {

                    # Assume they're using 1 & 2 as their allele "names"
                    if (!defined($LociAttributes{$Name}{Alleles}{$Allele})) {
                        $LociAttributes{$Name}{Alleles}{$Allele}{Order} = "$Allele";
                        push @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, $Allele;
                    }
                }

                # Translate the allele to a sequence number, this is required for bucket pattern matching
                $Allele = $LociAttributes{$Name}{Alleles}{$Allele}{Order};
            }

            # Pair-up the alleles
            if ($AlC & 1) {
                $Left = $Allele;
            } else {
                push @Pairs, $Left . " " . $Allele;
            }
        }

        if ($Type ne "BARE") {
            if (!$noparents) {
                $Pedigrees{$Ped}{$Ind}{Dad} = $Dad;
                $Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
            }
            if ($Type ne "PRE") {
                $Pedigrees{$Ped}{$Ind}{Kid1} = $Kid1;
                $Pedigrees{$Ped}{$Ind}{nPs}  = $nPs;    # Next paternal sibling
                $Pedigrees{$Ped}{$Ind}{nMs}  = $nMs;    # Next maternal sibling
                $Pedigrees{$Ped}{$Ind}{Prb}  = $Prb;    # Proband
            }
            $Pedigrees{$Ped}{$Ind}{Sex} = $Sex;
        } else {
            $Pedigrees{$Ped}{$Ind}{Sex} = 1;            # Bare case-control can be all male
        }
        $Pedigrees{$Ped}{$Ind}{Aff}    = $Aff;
        $Pedigrees{$Ped}{$Ind}{AlC}    = $AlC;          # Allele count (should always be the same)
        $Pedigrees{$Ped}{$Ind}{GtC}    = $GtC;          # Genotype count (how complete is genotyping)
        $Pedigrees{$Ped}{$Ind}{Mks}    = \@Pairs;
        $Pedigrees{$Ped}{$Ind}{OldFam} = $OldFam;
        $Pedigrees{$Ped}{$Ind}{OldInd} = $OldInd;
    }
    $PairCount = $AlC / 2;
    close IN;

    # Adopt single case/control individuals
    for my $Ped (sort keys %Pedigrees) {
        my $memberCount = scalar(keys %{ $Pedigrees{$Ped} });
        if ($memberCount == 1) {
            for my $Ind (keys %{ $Pedigrees{$Ped} }) {    # Yes, there's only one, but what individual ID?
                print "Adding parents for Pedigree $Ped, Individual $Ind\n";
                $Pedigrees{$Ped}{$Ind}{Prb} = 1;
                my $Dad = $Ind . "D";
                $Pedigrees{$Ped}{$Ind}{Dad} = $Dad;
                $Pedigrees{$Ped}{$Dad}{Sex} = 1;
                $Pedigrees{$Ped}{$Dad}{Dad} = $UnknownPerson;
                $Pedigrees{$Ped}{$Dad}{Mom} = $UnknownPerson;
                $Pedigrees{$Ped}{$Dad}{Aff} = $UnknownAffection;
                $Pedigrees{$Ped}{$Dad}{Prb} = 0;
                $Pedigrees{$Ped}{$Dad}{AlC} = $Pedigrees{$Ped}{$Dad}{GtC} = 0;
                $Pedigrees{$Ped}{$Dad}{LC}  = 1;
                $Pedigrees{$Ped}{$Dad}{Mks} = [ ("0 0") x ($Pedigrees{$Ped}{$Ind}{AlC} / 2) ];
                my $Mom = $Ind . "M";
                $Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
                $Pedigrees{$Ped}{$Mom}{Sex} = 2;
                $Pedigrees{$Ped}{$Mom}{Dad} = $UnknownPerson;
                $Pedigrees{$Ped}{$Mom}{Mom} = $UnknownPerson;
                $Pedigrees{$Ped}{$Mom}{Aff} = $UnknownAffection;
                $Pedigrees{$Ped}{$Mom}{Prb} = 0;
                $Pedigrees{$Ped}{$Mom}{AlC} = $Pedigrees{$Ped}{$Mom}{GtC} = 0;
                $Pedigrees{$Ped}{$Mom}{LC}  = 1;
                $Pedigrees{$Ped}{$Mom}{Mks} = [ ("0 0") x ($Pedigrees{$Ped}{$Ind}{AlC} / 2) ];
            }
        }
    }
    return $PairCount;
}

#####################################
# Derive allele frequencies from pedigree data when we don't have a config file.
#
sub deriveAlleleFrequencies {
    for my $i (0 .. $PairCount - 1) {
        $LociAttributes{ $Loci[ $i ] }{Type}         = "M";
        $LociAttributes{ $Loci[ $i ] }{Included}     = 1;
        my %HaploCounts = ();
        for my $Ped (keys %Pedigrees) {
            for my $Ind (keys %{ $Pedigrees{$Ped} }) {
                next if ($Pedigrees{$Ped}{$Ind}{Aff} != $Unaffected);
                my ($Left, $Right) = split /\s+/, $Pedigrees{$Ped}{$Ind}{Mks}[$i];
                $HaploCounts{$Left}++  if ($Left  ne AttributeMissing);
                $HaploCounts{$Right}++ if ($Right ne AttributeMissing);
            }
        }
#	print "HaploCounts for ".$Loci[ $i ]." are: ".Dumper(\%HaploCounts)."\n";
        my $PopSize = 0;
	for my $Allele (keys %HaploCounts) {
	    $PopSize += $HaploCounts{$Allele};
	}
	# Compute the ones for which we have occurrances
	for my $Allele (keys %HaploCounts) {
	    $LociAttributes{ $Loci[ $i ] }{Alleles}{$Allele}{Frequency} = $HaploCounts{$Allele} / $PopSize
		if ($PopSize != 0);
	}
	# Fill-in the rest with zero
	for my $Allele (@{ $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList} }) {
	    $LociAttributes{ $Loci[ $i ] }{Alleles}{$Allele}{Frequency} = '0'
		if (!defined($LociAttributes{ $Loci[ $i ] }{Alleles}{$Allele}{Frequency}));
	}
    }
}

#####################################
# Open and read the marker description companion file producing an ordered list of
# loci names (@Loci) and a hash for named loci attributes (%LociAttributes).
# Probably should fold names into the hash to get rid of @Loci.
#
# Note that Merlin treats 'A' as affectation status (binary), and 'M' as a
# quantitative trait. Since kelvin treats them the same, so do we.
#
# Added 'C' as liability class column. Return true if encountered.
#
sub loadCompanion {
    my $liability = 0;
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    my $Order  = 0;
    @Loci           = ();
    %LociAttributes = ();
    my @PedColUse = ();
    while (<IN>) {
        $LineNo++;
        print "Read $LineNo lines of $File\n" if (($LineNo % 1001) == 1000);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        my ($Type, $Name) = split /\s+/;
	push @PedColUse, $Type;
	die "Unknown locus type \"$Type\" at line $LineNo in marker description companion file $File\n"
	    if ($Type !~ /[MCAT]/);
	$liability = 1 if ($Type eq "C");
	if ($Type =~ /[M]/) {
	    push @Loci, $Name;
	    $LociAttributes{$Name}{Type}     = $Type;
	    $LociAttributes{$Name}{Included} = 1;
	}
    }
    close IN;
    return $liability;
}

#####################################
# Open and read the marker file to flesh-out the %LociAttributes hash. Named
# alleles get mapped to numbers so we can use our allele patterns when bucketizing.
#
sub loadFrequencies {
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo      = 0;
    my $Name        = "";
    my $AlleleCount = 0;
    while (<IN>) {
        $LineNo++;
        print "Read $LineNo lines of $File\n" if (($LineNo % 1001) == 1000);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        my @Tokens     = split /\s+/;
        my $RecordType = shift @Tokens;
        if ($RecordType eq "M") {
            $Name = shift @Tokens;
            die "Marker $Name at line $LineNo of $File not found in map file.\n" if (!defined($Map{$Name}{SAPos}));
            $AlleleCount = 0;
        } elsif ($RecordType eq "F") {    # List of unnamed allele frequencies
            # Dummy-up names so we know we'll always have them
            for (1 .. scalar(@Tokens)) {
                $AlleleCount++;
                $LociAttributes{$Name}{Alleles}{$AlleleCount}{Order} = $AlleleCount;
		push @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, $AlleleCount;
                $LociAttributes{$Name}{Alleles}{$AlleleCount}{Frequency} = shift(@Tokens);
            }
        } elsif ($RecordType eq "A") {    # Named allele and frequency
            $LociAttributes{$Name}{Alleles}{ $Tokens[0] }{Order} = ++$AlleleCount;
            push @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, $Tokens[0];
            $LociAttributes{$Name}{Alleles}{ $Tokens[0] }{Frequency} = $Tokens[1]
        } else {
            die "Unknown record type \"$RecordType\" at line $LineNo in marker file $File\n";
        }
    }
    close IN;
    return;
}

#####################################
# Add any implied biallelic marker alleles. This is pretty shakey as it only works with
# unnamed alleles.
#
sub addMissingAlleles {
    my $maf0 = 0; # True if we found missing minor allele frequencies
    for my $Name (@Loci) {
        next if ($LociAttributes{$Name}{Type}  =~ /^[AT]$/);
	if (defined($LociAttributes{$Name}{Alleles}{2}) && !defined($LociAttributes{$Name}{Alleles}{1})) {
	    $maf0 = 1;
	    $LociAttributes{$Name}{Alleles}{1}{Order} = 1;
	    unshift @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, '1';
	    $LociAttributes{$Name}{Alleles}{2}{Order} = 2;
	    $LociAttributes{$Name}{Alleles}{1}{Frequency} = '0';
	}
	if (defined($LociAttributes{$Name}{Alleles}{1}) && !defined($LociAttributes{$Name}{Alleles}{2})) {
	    $maf0 = 1;
	    $LociAttributes{$Name}{Alleles}{2}{Order} = 2;
	    push @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, '2';
	    $LociAttributes{$Name}{Alleles}{2}{Frequency} = '0';
	}
	if (!defined($LociAttributes{$Name}{Alleles}{1}) && !defined($LociAttributes{$Name}{Alleles}{2})) {
	    $maf0 = 1;
	    $LociAttributes{$Name}{Alleles}{1}{Order} = 1;
	    $LociAttributes{$Name}{Alleles}{2}{Order} = 2;
	    push @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, '1';
	    push @{ $LociAttributes{$Name}{Alleles}{OrderedList} }, '2';
	    $LociAttributes{$Name}{Alleles}{1}{Frequency} = '0';
	    $LociAttributes{$Name}{Alleles}{2}{Frequency} = '0';
	}	    
	@{ $LociAttributes{$Name}{Alleles}{OrderedList} } = sort @{ $LociAttributes{$Name}{Alleles}{OrderedList} };
    }
    return $maf0;
}

#####################################
# Open and read the map file for completeness' sake.
#
sub loadMap {
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    while (<IN>) {
        $LineNo++;
        print "Read $LineNo lines of $File\n" if (($LineNo % 1001) == 1000);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);          # Drop empty lines
        last if (/^\s*chr/i);    # Reached the header for map data
        ($MapFunction) = /^\s*mapfunction\s*=\s*(.*)$/i;
    }
    while (<IN>) {
        my ($Chr, $Marker, $SAPos, $MalePos, $FemalePos, $BPPos) = split /\s+/;
        $Map{$Marker}{Chr} = $Chr;
        $Map{$Marker}{SAPos} = $SAPos;
        $Map{$Marker}{MalePos} = $SAPos;
        $Map{$Marker}{FemalePos} = $SAPos;
        $Map{$Marker}{BPPos} = $SAPos;
    }
    close IN;
}

#####################################
# Check inter-file integrity, i.e. affectation and markers
#
sub checkIntegrity {
    for my $Ped (sort keys %Pedigrees) {
        my $UnknownAffectionCount = 0;
        my $UnaffectedCount       = 0;
        my $AffectedCount         = 0;
        my $MagicCounts           = 0;
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            $UnknownAffectionCount++ if ($Pedigrees{$Ped}{$Ind}{Aff} == $UnknownAffection);
            $UnaffectedCount++       if ($Pedigrees{$Ped}{$Ind}{Aff} == $Unaffected);
            $AffectedCount++         if ($Pedigrees{$Ped}{$Ind}{Aff} == $Affected);
            $MagicCounts++           if ($Pedigrees{$Ped}{$Ind}{Aff} =~ /(88\.88|99\.99|NaN|Inf)/i);
            my @Pairs = @{ $Pedigrees{$Ped}{$Ind}{Mks} };
            for my $i (0 .. $PairCount - 1) {
                my ($Left, $Right) = split /\s/, $Pairs[$i];
                die "Pedigree $Ped, individual $Ind Marker $i (" . $Loci[ $i ] . ") allele $Left too large.\n"
                  if ($Left > scalar( @{ $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList} }));
                die "Pedigree $Ped, individual $Ind Marker $i (" . $Loci[ $i ] . ") allele $Right too large.\n"
                  if ($Right > scalar( @{ $LociAttributes{ $Loci[ $i ] }{Alleles}{OrderedList} }));
                if ((defined($Directives{XC}) || $XC)) {
                    die "Pedigree $Ped, male $Ind is not homozygous for marker " . $Loci[ $i ] . "\n"
                      if (($Pedigrees{$Ped}{$Ind}{Sex} == 1) && ($Left != $Right));
                }
            }
        }
        if (defined($Directives{QT})) {
            if ($UnknownAffectionCount + $UnaffectedCount + $AffectedCount < $MagicCounts) {
                print "Warning! Your QT analysis for pedigree $Ped has Unk/UnA/Aff of "
                  . "$UnknownAffectionCount/$UnaffectedCount/$AffectedCount\nout of "
                  . scalar(keys %{ $Pedigrees{$Ped} })
                  . " individuals and "
                  . $MagicCounts
                  . " default value(s) (any of 88.88/99.99/NaN/Inf).\n";
            }
        } else {

            # Must be DT
            die "Your DT analysis for pedigree $Ped has Unk/UnA/Aff of "
              . "$UnknownAffectionCount/$UnaffectedCount/$AffectedCount out of "
              . scalar(keys %{ $Pedigrees{$Ped} })
              . " individuals\n"
              if ($UnknownAffectionCount + $UnaffectedCount + $AffectedCount != scalar(keys %{ $Pedigrees{$Ped} }));
        }
    }

}

#####################################
# Check relations and genders (parent, siblings present, not self-parenting, correct genders)
#
sub checkRelations {

    my $Type = shift();
    for my $Ped (sort keys %Pedigrees) {
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            if ($Dad ne $UnknownPerson) {
                die "In pedigree $Ped, $Ind is own father!\n" if ($Dad eq $Ind);
                die "In pedigree $Ped, Father $Dad not found for individual $Ind!\n"
                  if (!defined($Pedigrees{$Ped}{$Dad}));
                die "In pedigree $Ped, Father $Dad is not male!\n" if ($Pedigrees{$Ped}{$Dad}{Sex} != 1);
            }
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if ($Mom ne $UnknownPerson) {
                die "In pedigree $Ped, $Ind is own mother!\n" if ($Mom eq $Ind);
                die "In pedigree $Ped, Mother $Mom not found for individual $Ind!\n"
                  if (!defined($Pedigrees{$Ped}{$Mom}));
                die "In pedigree $Ped, Mother $Mom is not female!\n" if ($Pedigrees{$Ped}{$Mom}{Sex} != 2);
            }
            if ($Type eq "POST") {
                my $Kid1 = $Pedigrees{$Ped}{$Ind}{Kid1};
                if ($Kid1 ne $UnknownPerson) {
                    die "In pedigree $Ped, $Ind is own sibling!\n" if ($Kid1 eq $Ind);
                    die "In pedigree $Ped, first child $Kid1 missing for individual $Ind!\n"
                      if (!defined($Pedigrees{$Ped}{$Kid1}));
                }
            }
        }
    }
}

sub numerically { $a <=> $b }

sub numericIsh {
    if (($a . $b) =~ /^\d+$/) {
        $a <=> $b;
    } else {
        $a cmp $b;
    }
}

#####################################
# Shamelessly stolen from an example on the Internet at http://snippets.dzone.com/posts/show/99
sub expand {
    my $Range = shift;
    my @Result;
    $Range =~ s/[^\d\-\,]//gs;    #remove extraneous characters
    my @Items = split(/,/, $Range);
    foreach (@Items) {
        m/^\d+$/ and push(@Result, $_) and next;
        my ($Start, $Finish) = split /-/;
        push(@Result, ($Start .. $Finish)) if $Start < $Finish;
    }
    return @Result;
}

#####################################
# Start discovering statistics that might affect performance.
#
sub perfStats {

    my %CountsLabels = (
        MsgMkr => "Missing Markers",
        MsgAff => "Missing Affected Status",
        MsgSex => "Missing Sex",
        Fdrs   => "Founders",
        NoFdrs => "Non-Founders"
    );
    my %PedCat = ();    # Hash of counts for labelled categories indexed by pedigree
    my @PedSiz = ();    # List of counts of individuals in pedigree
    my @PedNuk = ();    # List of counts of nuclear families in pedigrees
    my $Avg;
    for my $Ped (keys %Pedigrees) {
        for my $Label (keys %CountsLabels) {
            $PedCat{$Ped}{$Label} = 0;
        }
        my %Seen = ();    # Parental pairs seen (i.e. nuclear families seen)
        for my $Ind (sort keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom ne $UnknownPerson) && ($Dad ne $UnknownPerson)) {
                $Seen{ $Mom . "+" . $Dad }++;
            }
            if (($Dad eq $UnknownPerson) && ($Mom eq $UnknownPerson)) {
                $PedCat{$Ped}{Fdrs}++;
            } else {
                $PedCat{$Ped}{NoFdrs}++;
            }
            if ($Pedigrees{$Ped}{$Ind}{GtC} == 0) { $PedCat{$Ped}{MsgMkr}++; }
            if ($Pedigrees{$Ped}{$Ind}{Aff} == $UnknownAffection) { $PedCat{$Ped}{MsgAff}++; }
            if ($Pedigrees{$Ped}{$Ind}{Sex} ne AttributeMissing) { $PedCat{$Ped}{MsgSex}++; }
        }
        push @PedSiz, scalar(keys(%{ $Pedigrees{$Ped} }));
        push @PedNuk, scalar(keys(%Seen));
    }

    @PedSiz = sort numerically @PedSiz;
    print "Number of pedigrees:" . @PedSiz . "\n";
    $Avg = sprintf("%.2f", (sum @PedSiz) / @PedSiz);
    @PedNuk = sort numerically @PedNuk;
    print "...Sizes: min:" . $PedSiz[0] . " max:" . $PedSiz[-1] . " med:" . $PedSiz[ $#PedSiz / 2 ] . " avg: $Avg\n";
    $Avg = sprintf("%.2f", (sum @PedNuk) / @PedNuk);
    print "...Nuclear family counts: min:"
      . $PedNuk[0] . " max:"
      . $PedNuk[-1] . " med:"
      . $PedNuk[ $#PedNuk / 2 ]
      . " avg: $Avg\n";

    for my $Label (keys %CountsLabels) {
        my @Counts = ();
        for my $Ped (keys %Pedigrees) {
            push @Counts, $PedCat{$Ped}{$Label};
        }
        @Counts = sort numerically @Counts;
        print "..." . $CountsLabels{$Label} . ": ";
        $Avg = sprintf("%.2f", (sum @Counts) / @Counts);
        print "min:" . $Counts[0] . " max:" . $Counts[-1] . " med:" . $Counts[ $#Counts / 2 ] . " avg: $Avg\n";
    }

    # Depth statistics -- the number of direct ancestors
    for my $Ped (keys %Pedigrees) {
        @Ancestors = ();
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            push @Ancestors, &countAncestors($Ped, $Ind, 0);
        }
        @Ancestors = sort numerically @Ancestors;
        $Avg = sprintf("%.2f", (sum @Ancestors) / @Ancestors);
        print "Pedigree $Ped Depth: max:" . $Ancestors[-1] . " med:" . $Ancestors[ $#Ancestors / 2 ] . " avg: $Avg\n";
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
      if (defined($Directives{DA}) && ($Directives{DA}[0] != 2));

    if (defined($Directives{DK})) {

        # Integration (DK) analysis limitations
    } else {

        # Iterative (classic Kelvin) analysis limitations
    }
}

1
