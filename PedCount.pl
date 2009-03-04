#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(sum);
use Data::Dumper;

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

# Command line option flags
my $config = 0; my $pre = 0; my $post = 0; my $noparents = 0; my $XC = 0; 
my $bare = 0; my $count = 0; my $write = "unspecified"; my $loops = 0; my $stats = 0;
my $split = 0; my $nokelvin = 0; my @include = (); my @exclude = ();
my $WritePrefix = "PC";

# Permanent defaults
use constant AttributeMissing => "0";    # For marker alleles and Sex

# Defaults to be overridden by configuration file directives
my $pedFile          = "pedpost.dat";
my $mapFile          = "mapfile.dat";
my $companionFile    = "datafile.dat";
my $markersFile      = "markers.dat";    # Default input files
my $UnknownAffection = 0;
my $Unaffected       = 1;
my $Affected         = 2;                # Default affection indicators
my $UnknownPerson    = "0";              # Default unknown person indicator
my $MapFunction = "Kosambi, bless his heart!";

# Read/calculated results
my %Pedigrees;                           # Pedigrees as loaded
my %Directives;                          # Directives as loaded
my $PairCount = 0;                       # Last pedigree count of marker pairs
my @Loci;                                # Ordered loci name list from companion file
my %LociAttributes;                      # Loci attributes from companion and marker files
my %Map; # Loci on the map and other map attributes

# Nuisances to fix
my @Depths;                              # Referenced recursively and I'm fuddled
my @Ancestors;                           # ditto
my $ShortestLoop = "";                   # Just batted around too much to monkey with right now.

# Known directives as well as dispatch routine if needed. I could avoid the NoActions, but
# I prefer to be explicit.
my %KnownDirectives = (
    AL    => \&NoAction,
    AM    => \&NoAction,
    AS    => \&dirAS,
    CC    => \&dirCF,
    CF    => \&dirCF,
    DA    => \&NoAction,
    DD    => \&NoAction,
    DF    => \&NoAction,
    DK    => \&NoAction,
    DT    => \&NoAction,
    Dd    => \&NoAction,
    GF    => \&NoAction,
    HE    => \&NoAction,
    LC    => \&NoAction,
    LD    => \&NoAction,
    LOG   => \&NoAction,
    MK    => \&NoAction,
    MM    => \&NoAction,
    MP    => \&NoAction,
    MX    => \&NoAction,
    P1    => \&NoAction,
    PD    => \&NoAction,
    PE    => \&NoAction,
    PF    => \&NoAction,
    QT    => \&dirQT,
    SA    => \&NoAction,
    SS    => \&NoAction,
    TL    => \&NoAction,
    TM    => \&NoAction,
    TP    => \&NoAction,
    TT    => \&NoAction,
    Th    => \&NoAction,
    T_MIN => \&NoAction,
    T_MAX => \&NoAction,
    UP    => \&dirUP,
    XC    => \&dirXC,
    dD    => \&NoAction,
    dd    => \&NoAction,
);

#####################################
#
sub dirAS {
    $UnknownAffection = $Directives{AS}[0];
    $Unaffected       = $Directives{AS}[1];
    $Affected         = $Directives{AS}[2];
}

#####################################
#
sub dirCF {
    die "Cannot generate counts for configuration that already has them."
      if ($count);
}

#####################################
#
sub dirQT {

    # If it's not already explicitly specified, set QT default for AS
    if (!defined($Directives{AS})) {
        $UnknownAffection = -99.99;
        $Unaffected       = -88.88;
        $Affected         = 88.88;
    }
}

#####################################
#
sub dirUP {
    $UnknownPerson = $Directives{UP}[0];
}

#####################################
#
sub dirXC {
    print "meh.\n";
}

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
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open $File\n";
    my $LineNo = 0;
    %Directives = ();
    while (<IN>) {
        $LineNo++;
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        my @Parameters = split;
        my $Directive  = shift @Parameters;
        die "Configuration line $LineNo: \"$Directive\" is not a known configuration directive."
          if (!defined($KnownDirectives{$Directive}));
        $Directives{$Directive} = \@Parameters;
        &{ $KnownDirectives{$Directive} };
    }
}

#####################################
# Open and assess a pedigree file to determine processing required. Return a type, one of
# POST, PRE or BARE.
#
sub assessPedigree {
    my $File = shift();
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
        my $MarkerColumns = scalar(@Loci) * 2;
        if ($pre) {
            $TotalColumns += 2 if ($noparents);    # Make up for no parents here
                  # pre-MAKEPED is Ped Ind Dad Mom Sex Aff Pairs, i.e. some even number > 6
            die "Invalid number of columns in pre-MAKEPED pedigree file $File."
              if (($TotalColumns < 8) || ($TotalColumns & 1));
            die
              "Inconsistent column counts between companion datafile $companionFile and pre-MAKEPED pedigree file $pedFile."
              if ($config && (($TotalColumns - $MarkerColumns) != 6));
            return "PRE";
        } elsif ($post) {
            $TotalColumns += 2 if ($noparents);    # Make up for no parents here
                  # post-MAKEPED is Ped Ind Dad Mom Kid1 nPs nMs Sex Prb Aff Pairs, i.e. some even number > 10
            die "Invalid number of columns in post-MAKEPED pedigree file $File."
              if (($TotalColumns < 12) || ($TotalColumns & 1));
            die
              "Inconsistent column counts between companion datafile $companionFile and post-MAKEPED pedigree file $pedFile."
              if ($config && (($TotalColumns - $MarkerColumns) != 10));
            return "POST";
        } elsif ($bare) {

            # Bare individuals for CC are $Ped/$Ind Aff Pairs, i.e. some even number > 2
            die "Invalid number of columns in bare pedigree file $File."
              if (($TotalColumns < 4) || ($TotalColumns & 1));
            die
              "Inconsistent column counts between companion datafile $companionFile and post-MAKEPED pedigree file $pedFile."
              if ($config && (($TotalColumns - $MarkerColumns) != 2));
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
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    %Pedigrees = ();
    my $AlC = 0; # Number of alleles
    my $GtC = 0; # Number of genotyped alleles
    my ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Alleles);

    while (<IN>) {
        $LineNo++;
	print "At line $LineNo of $File\n" if (($LineNo % 1024) == 1023);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        if ($Type eq "POST") {
            ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Alleles) = split /\s+/;
            die "Pedigree line $LineNo: proband must be numeric, not \"$Prb\"." if (!($Prb =~ /[0-9]/));
        } elsif ($Type eq "PRE") {
            ($Ped, $Ind, $Dad, $Mom, $Sex, $Aff, @Alleles) = split /\s+/;
        } elsif ($Type eq "BARE") {
            ($Ind, $Aff, @Alleles) = split /\s+/;
            $Ped = $Ind;
        } else {
            die "Unhandled pedigree type \"$Type\".";
        }

        # Validate everything we've got so far
        die "Pedigree line $LineNo: sex must be 1 or 2, not \"$Sex\"." if ((!$bare) && !($Sex =~ /[12]/));
        die
          "Pedigree line $LineNo: affection status must be $UnknownAffection, $Unaffected, or $Affected, not \"$Aff\"."
          if ( (!defined($Directives{QT}))
            && (($Aff != $UnknownAffection) && ($Aff != $Unaffected) && ($Aff != $Affected)));

	# Handle marker pairs
        my ($OldFam, $OldInd, $Left);
        $AlC = 0;
        $GtC = 0;
        my @Pairs = ();
        for my $Allele (@Alleles) {
            if ($Allele eq "Ped:") {
                $OldFam = $Alleles[ $AlC + 1 ];
                $OldInd = $Alleles[ $AlC + 3 ];
                last;
            }
            $AlC++;
	    if ($Allele ne AttributeMissing) {
		$GtC++; # Keep track of how much genotypic information we have for this individual
		my $Name = "";
		if (!$config) {
		    # We need to fake-up @Loci and %LociAttributes so we can still translate from 
		    # named to sequenced alleles. After all, alleles could be "2" and "foo".
		    $Name = int(($AlC + 1.5) / 2); # Offset == Name
		    if (!defined($LociAttributes{$Name}{Alleles}{$Allele})) {
			# If this allele is not already know, it's the next one sequentially
			$LociAttributes{$Name}{Alleles}{$Allele} = scalar(keys %{ $LociAttributes{$Name}{Alleles} }) + 1;
		    }
		} else {
		    $Name = $Loci[int(($AlC + 1.5) / 2)]; # Integer division, thank you
		}
		# Translate the allele to a sequence number
		$Allele = $LociAttributes{$Name}{Alleles}{$Allele};
	    }

	    # Pair-up the alleles
            if ($AlC & 1) {
                $Left = $Allele;
            } else {
                push @Pairs, $Left . " " . $Allele;
            }
        }

        if (!$bare) {
            if (!$noparents) {
                $Pedigrees{$Ped}{$Ind}{Dad} = $Dad;
                $Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
            }
            if (!$pre) {
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
                $Pedigrees{$Ped}{$Dad}{Mks} = [ ("0 0") x ($Pedigrees{$Ped}{$Ind}{AlC} / 2) ];
                my $Mom = $Ind . "M";
                $Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
                $Pedigrees{$Ped}{$Mom}{Sex} = 2;
                $Pedigrees{$Ped}{$Mom}{Dad} = $UnknownPerson;
                $Pedigrees{$Ped}{$Mom}{Mom} = $UnknownPerson;
                $Pedigrees{$Ped}{$Mom}{Aff} = $UnknownAffection;
                $Pedigrees{$Ped}{$Mom}{Prb} = 0;
                $Pedigrees{$Ped}{$Mom}{AlC} = $Pedigrees{$Ped}{$Mom}{GtC} = 0;
                $Pedigrees{$Ped}{$Mom}{Mks} = [ ("0 0") x ($Pedigrees{$Ped}{$Ind}{AlC} / 2) ];
            }
        }
    }
}

#####################################
# Derive a list of loci and marker alleles and frequencies from the pedigree itself
# and store the names in the global list @Loci and attributes in the global two-
# dimensional hash %LociAttributes.
#
sub deriveLociAndAttributes {
    @Loci                        = ();
    %LociAttributes              = ();
    $Loci[0]                     = "Trait";
    $LociAttributes{Trait}{Type} = "T";
    $LociAttributes{Trait}{Included} = 0;
    for my $i (0 .. $PairCount - 1) {
        $Loci[ $i + 1 ] = sprintf("M%04d", $i + 1);
        $LociAttributes{ $Loci[ $i + 1 ] }{Type} = "M";
        $LociAttributes{ $Loci[ $i + 1 ] }{Included} = 1;
        my @HaploCounts = (0); # List works because they must be numeric (relative position of frequency in marker file)
        for my $Ped (keys %Pedigrees) {
            for my $Ind (keys %{ $Pedigrees{$Ped} }) {
                next if ($Pedigrees{$Ped}{$Ind}{Aff} != $Unaffected);
                my ($Left, $Right) = split /\s+/, $Pedigrees{$Ped}{$Ind}{Mks}[$i];
                $HaploCounts[$Left]++  if ($Left ne AttributeMissing);
                $HaploCounts[$Right]++ if ($Right ne AttributeMissing);
            }
        }
        my $PopSize = sum @HaploCounts;
        my @Tokens  = ();
        for $i (1 .. scalar(@HaploCounts) - 1) {
            push @Tokens, $HaploCounts[$i] / $PopSize;
        }
        $LociAttributes{ $Loci[ $i + 1 ] }{Frequencies} = [@Tokens];
    }
}

#####################################
# Open and read the marker description companion file producing an ordered list of
# loci names (@Loci) and a hash for named loci attributes (%LociAttributes).
# Probably should fold names into the hash to get rid of @Loci.
#
sub loadCompanion {
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0; my $Order = 0;
    @Loci           = ();
    %LociAttributes = ();
    while (<IN>) {
        $LineNo++;
	print "At line $LineNo of $File\n" if (($LineNo % 1024) == 1023);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        my ($Type, $Name) = split /\s+/;
        die "Unknown locus type \"$Type\" at line $LineNo in marker description companion file $File\n"
          if (($Type ne "T") and ($Type ne "A") and ($Type ne "M"));
        push @Loci, $Name;
        $LociAttributes{$Name}{Type} = $Type;
        $LociAttributes{$Name}{Included} = 1;
    }
    close IN;
}

#####################################
# Open and read the marker file to flesh-out the %LociAttributes hash. Named
# alleles get mapped to numbers so we can use our allele patterns when bucketizing.
#
sub loadMarkers {
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    my $Name   = "";
    my $AlleleCount = 0;
    while (<IN>) {
        $LineNo++;
	print "At line $LineNo of $File\n" if (($LineNo % 1024) == 1023);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        my @Tokens     = split /\s+/;
        my $RecordType = shift @Tokens;
        if ($RecordType eq "M") {
            $Name = shift @Tokens;
	    die "Marker $Name at line $LineNo of $File not found in map file.\n" if (!defined($Map{$Name}{Pos}));
	    $AlleleCount = 0;
        } elsif ($RecordType eq "F") { # List of unnamed allele frequencies
            if (defined($LociAttributes{$Name}{Frequencies})) {
                $LociAttributes{$Name}{Frequencies} = [ (@{ $LociAttributes{$Name}{Frequencies} }, @Tokens) ];
            } else {
                $LociAttributes{$Name}{Frequencies} = [@Tokens];
            }
	    # Dummy-up names so we know we'll always have them
	    for (1..scalar(@Tokens)) {
		$AlleleCount++;
		$LociAttributes{$Name}{Alleles}{$AlleleCount} = $AlleleCount;
	    }
        } elsif ($RecordType eq "A") { # Named allele and frequency
	    push @{ $LociAttributes{$Name}{Frequencies} }, $Tokens[1];
	    $LociAttributes{$Name}{Alleles}{$Tokens[0]} = ++$AlleleCount;
	} else {
            die "Unknown record type \"$RecordType\" at line $LineNo in marker file $File\n";
        }
    }
    close IN;
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
	print "At line $LineNo of $File\n" if (($LineNo % 1024) == 1023);
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
	last if (/^\s*chr/i); # Reached the header for map data
	($MapFunction) = /^\s*mapfunction\s*=\s*(.*)$/i;
    }
    while (<IN>) {
	my ($Chr, $Marker, $Pos) = split /\s+/;
	$Map{$Marker}{Chr} = $Chr;
	$Map{$Marker}{Pos} = $Pos;
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
                die "Pedigree $Ped, individual $Ind Marker $i (" . $Loci[ $i + 1 ] . ") allele $Left too large.\n"
                  if ($Left > scalar(@{ $LociAttributes{ $Loci[ $i + 1 ] }{Frequencies} }));
                die "Pedigree $Ped, individual $Ind Marker $i (" . $Loci[ $i + 1 ] . ") allele $Right too large.\n"
                  if ($Right > scalar(@{ $LociAttributes{ $Loci[ $i + 1 ] }{Frequencies} }));
                if ((defined($Directives{XC}) || $XC)) {
                    die "Pedigree $Ped, male $Ind is not homozygous for marker " . $Loci[ $i + 1 ] . "\n"
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
    $Range =~ s/[^\d\-\,]//gs; #remove extraneous characters
    my @Items = split(/,/,$Range);    
    foreach (@Items){                 
	m/^\d+$/ and push(@Result,$_) and next;
	my ($Start, $Finish) = split /-/;
	push(@Result,($Start .. $Finish)) if $Start < $Finish;    }                                 
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
            if ($Pedigrees{$Ped}{$Ind}{GtC} == 0)                 { $PedCat{$Ped}{MsgMkr}++; }
            if ($Pedigrees{$Ped}{$Ind}{Aff} == $UnknownAffection) { $PedCat{$Ped}{MsgAff}++; }
            if ($Pedigrees{$Ped}{$Ind}{Sex} ne AttributeMissing)  { $PedCat{$Ped}{MsgSex}++; }
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
#
# Adapted from a trios-only version by John Burian.
#
# Leverages pedigree inheritance patterns to reduce computational
# complexity. Multiple pedigrees that conform to a single inheritance pattern
# are converted to a single pedigree and a weight (count) for each marker.
# Currently works only for two-point analysis because it would be a real
# trick for there to be many pedigrees where multiple markers all fit the
# same inheritance pattern. 
#
# The traditional analyses for which this was designed are:
# 1. Case/Control, where unrelated genotyped and phenotyped individuals are 
# divided into cases and controls and fall into one of 6 inheritance pattern
# buckets normally, 10 for X chromosome analysis.
# 2. Trios are 3-member pedigrees with one affected child both parents treated
# as if unaffected. They fall into one of 30 inheritance pattern buckets
# normally, and 44 for X chromosome analysis.
# 3. ASPs (Affected Sibling Pairs) are essentially expanded trios, i.e. 4-member 
# pedigrees with two affected children and both parents in 
# any affectation category, which will fall first into one of the 3-member 
# trio buckets for the first child, and then a second 3-member bucket
# with the same parentage for the second child, e.g. T30 and then T18. This 
# is used to identify a unique 4-person bucket by ordering the 3-member
# numbers, e.g. T18-T30, because a family falling into T30 and then T18 is the
# same as one falling into T18 and then T30. This extends to any number
# of children of a single-generation family.
#
# Affectation status is handled by concatenating it to either the parental or 
# child portion of the bucket name. Gender is handled for XC
# analysis by the same approach, i.e. concatenating it to either the parental
# or child portion of the bucket name. Ultimately we end-up with bucket names
# generated from the genotypic, phenotypic (and maybe even gender) attributes
# of the entire nuclear family. While enumerating every possible bucket for all
# possible nuclear family combinations would be exhaustive, we use Perl hashes
# to produce only the buckets needed. When we create a bucket, we keep track of
# the last pedigree that fit into it so we can use that pedigree's genotypic and
# phenotypic information to create the template pedigree for the bucket. We could
# decode the bucket name and achieve the same goal, but this is easier and more
# reliable.
#
sub bucketizePedigrees {

    # When determining non-XC buckets, alleles count but phase and parent doesn't (i.e.
    # 11+11=11 != 22+22=22, but 11+12=12 == 11+21=12 and 11+12=12 == 12+11=12).

    my %TrioBuckets = (    # Mom, then Dad, then the child, but it doesn't matter
        '0 0' => {
            '0 0' => {
                '0 0' => 'T30',    # 30  0 0  0 0  0 0
                '1 1' => 'T18',    # 18  0 0  0 0  1 1 (case11, control11)
                '1 2' => 'T19',    # 19  0 0  0 0  1 2 (case12, control12)
                '2 2' => 'T20',    # 20  0 0  0 0  2 2 (case22, control22)
            },
            '1 1' => {
                '0 0' => 'T24',    # 24  0 0  1 1  0 0
                '1 1' => 'T05',    #  5  0 0  1 1  1 1
                '1 2' => 'T06',    #  6  0 0  1 1  1 2
            },
            '1 2' => {
                '0 0' => 'T27',    # 27  0 0  1 2  0 0
                '1 1' => 'T12',    # 12  0 0  1 2  1 1
                '1 2' => 'T13',    # 13  0 0  1 2  1 2
                '2 2' => 'T14',    # 14  0 0  1 2  2 2
            },
            '2 2' => {
                '0 0' => 'T29',    # 29  0 0  2 2  0 0
                '1 2' => 'T16',    # 16  0 0  2 2  1 2
                '2 2' => 'T17',    # 17  0 0  2 2  2 2
            },
        },
        '1 1' => {
            '0 0' => {
                '0 0' => 'T24',    # 24  1 1  0 0  0 0
                '1 1' => 'T05',    #  5  1 1  0 0  1 1
                '1 2' => 'T06',    #  6  1 1  0 0  1 2
            },
            '1 1' => {
                '0 0' => 'T21',    # 21  1 1  1 1  0 0
                '1 1' => 'T01',    #  1  1 1  1 1  1 1
            },
            '1 2' => {
                '0 0' => 'T22',    # 22  1 1  1 2  0 0
                '1 1' => 'T02',    #  2  1 1  1 2  1 1
                '1 2' => 'T03',    #  3  1 1  1 2  1 2
            },
            '2 2' => {
                '0 0' => 'T23',    # 23  1 1  2 2  0 0
                '1 2' => 'T04',    #  4  1 1  2 2  1 2
            },
        },
        '1 2' => {
            '0 0' => {
                '0 0' => 'T27',    # 27  1 2  0 0  0 0
                '1 1' => 'T12',    # 12  1 2  0 0  1 1
                '1 2' => 'T13',    # 13  1 2  0 0  1 2
                '2 2' => 'T14',    # 14  1 2  0 0  2 2
            },
            '1 1' => {
                '0 0' => 'T22',    # 22  1 2  1 1  0 0
                '1 1' => 'T02',    #  2  1 2  1 1  1 1
                '1 2' => 'T03',    #  3  1 2  1 1  1 2
            },
            '1 2' => {
                '0 0' => 'T25',    # 25  1 2  1 2  0 0
                '1 1' => 'T07',    #  7  1 2  1 2  1 1
                '1 2' => 'T08',    #  8  1 2  1 2  1 2
                '2 2' => 'T09',    #  9  1 2  1 2  2 2
            },
            '2 2' => {
                '0 0' => 'T26',    # 26  1 2  2 2  0 0
                '1 2' => 'T10',    # 10  1 2  2 2  1 2
                '2 2' => 'T11',    # 11  1 2  2 2  2 2
            },
        },
        '2 2' => {
            '0 0' => {
                '0 0' => 'T29',    # 29  2 2  0 0  0 0
                '1 2' => 'T16',    # 16  2 2  0 0  1 2
                '2 2' => 'T17',    # 17  2 2  0 0  2 2
            },
            '1 1' => {
                '0 0' => 'T23',    # 23  2 2  1 1  0 0
                '1 2' => 'T04',    #  4  2 2  1 1  1 2
            },
            '1 2' => {
                '0 0' => 'T26',    # 26  2 2  1 2  0 0
                '1 2' => 'T10',    # 10  2 2  1 2  1 2
                '2 2' => 'T11',    # 11  2 2  1 2  2 2
            },
            '2 2' => {
                '0 0' => 'T28',    # 28  2 2  2 2  0 0
                '2 2' => 'T15',    # 15  2 2  2 2  2 2
            },
        },
    );

    # This bucket hash is a bit misleading because it uses the full-genotype 
    # notation even though the individual might be male and therefore have 
    # only one allele for the X chromosome.
    # When determining XC buckets, alleles and parent counts but phase doesn't (i.e.
    # 11+22=22 != 22+11=22 and 11+11=11 != 22+22=22, but 12+11=12 == 21+11=12).

    my %XCTrioBuckets = (          # Mom, then Dad, then the child, and it matters!
        '0 0' => {
            '0 0' => {
                '0 0' => 'X01',
                '1 1' => 'X02',    # XC case/control male or female (case11, control11, case1, control1)
                '1 2' => 'X03',    # XC case/control female-only (case12, control12)
                '2 2' => 'X04',    # XC case/control male or female (case22, control22, case2, control2)
            },
            '1 1' => {
                '0 0' => 'X05',
                '1 1' => 'X06',
                '1 2' => 'X07',    # Female
                '2 2' => 'X08',
            },
            '2 2' => {
                '0 0' => 'X09',
                '1 1' => 'X10',
                '1 2' => 'X11',    # Female
                '2 2' => 'X12',
            },
        },
        '1 1' => {
            '0 0' => {
                '0 0' => 'X13',
                '1 1' => 'X14',
                '1 2' => 'X15',    # Female
                '2 2' => 'X16',
            },
            '1 1' => {
                '0 0' => 'X17',
                '1 1' => 'X18',
            },
            '2 2' => {
                '0 0' => 'X18',
                '1 1' => 'X20',
                '1 2' => 'X21',    # Female
                '2 2' => 'X22',
            },
        },
        '1 2' => {
            '0 0' => {
                '0 0' => 'X23',
                '1 1' => 'X24',
                '1 2' => 'X25',    # Female
                '2 2' => 'X26',
            },
            '1 1' => {
                '0 0' => 'X27',
                '1 1' => 'X28',
                '1 2' => 'X29',    # Female
                '2 2' => 'X30',
            },
            '2 2' => {
                '0 0' => 'X31',
                '1 1' => 'X32',
                '1 2' => 'X33',    # Female
                '2 2' => 'X34',
            },
        },
        '2 2' => {
            '0 0' => {
                '0 0' => 'X35',
                '1 1' => 'X36',
                '1 2' => 'X37',    # Female
                '2 2' => 'X38',
            },
            '1 1' => {
                '0 0' => 'X39',
                '1 1' => 'X40',
                '1 2' => 'X41',    # Female
                '2 2' => 'X42',
            },
            '2 2' => {
                '0 0' => 'X43',
                '2 2' => 'X44',
            },
        },
    );

    my $Type = shift();            # Pedigree type for writing
    my $Prefix = shift(); # Uniqifying (what a word!) prefix for files

    # Verify that this is a 2pt analysis (default, so look for multipoint directives)
    die "Generation of counts not permitted for a multipoint analysis.\n"
      if (defined($Directives{SA}) || defined($Directives{SS}));

    # Verify that all markers are present and only biallelic...
    for my $Name (@Loci) {
	next if (!$LociAttributes{$Name}{Included});
        next if ($LociAttributes{$Name}{Type} eq "T");
        next if ($LociAttributes{$Name}{Type} eq "A");
        die "No allele information found for marker $Name for count generation.\n"
          if (!defined($LociAttributes{$Name}{Frequencies}));
        die "Marker $Name not biallelic, not permitted for count generation.\n"
          if (scalar(@{ $LociAttributes{$Name}{Frequencies} }) != 2);
    }

    my %Buckets   = ();    # Fully-funkified bucket names with encoded everything
    my @Skippies  = ();    # Pedigrees copied on thru without bucketization (affects stats)
    my %Templates = ();    # Template pedigrees
    my $PedSeq    = 1;     # Template pedigree ID to keep them short

    # Look at each family...
    for my $Ped (sort numericIsh keys %Pedigrees) {

        my $memberCount = scalar(keys %{ $Pedigrees{$Ped} });

        # Qualify the family for inclusion in trio buckets by
        # verifying depth of 1 while building a parental
        # affectation prefix so we can do more
        # than expected (i.e. handle any nuclear families)
        my $PAP = "";
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Dad eq $UnknownPerson) && ($Mom eq $UnknownPerson)) {
                if ($Pedigrees{$Ped}{$Ind}{Sex} == 1) {
                    $PAP = $Pedigrees{$Ped}{$Ind}{Aff} . $PAP;
                } else {
                    $PAP = $PAP . $Pedigrees{$Ped}{$Ind}{Aff};
                }
            } else {
                if (   ($Pedigrees{$Ped}{$Dad}{Dad} ne $UnknownPerson)
                    || ($Pedigrees{$Ped}{$Dad}{Mom} ne $UnknownPerson)
                    || ($Pedigrees{$Ped}{$Mom}{Dad} ne $UnknownPerson)
                    || ($Pedigrees{$Ped}{$Mom}{Mom} ne $UnknownPerson)) {
                    $PAP = "";
                    last;
                }
            }
        }
        if ($PAP eq "") {
            print "Will copy multi-generation pedigree $Ped intact.\n";
            push @Skippies, $Ped;
            next;
        }

        # Get the family genotype bucket for each marker pair
        for my $i (0 .. $PairCount - 1) {
            my @bucketList = ();

            # Get a trio bucket for each child in the family
            for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {

                # Skip parents
                my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
                if ($Dad ne $UnknownPerson) {
                    my $DadAlleles = $Pedigrees{$Ped}{$Dad}{Mks}[$i];
                    ($DadAlleles eq '2 1') and $DadAlleles = '1 2';
                    my $ChildAlleles = $Pedigrees{$Ped}{$Ind}{Mks}[$i];
                    ($ChildAlleles eq '2 1') and $ChildAlleles = '1 2';

                    my $Mom        = $Pedigrees{$Ped}{$Ind}{Mom};
                    my $MomAlleles = $Pedigrees{$Ped}{$Mom}{Mks}[$i];
                    ($MomAlleles eq '2 1') and $MomAlleles = '1 2';

#                        print "Get trio bucket for $MomAlleles $DadAlleles $ChildAlleles\n";
                    my $TrioBucket;
                    if (defined($Directives{XC}) || $XC) {
                        $TrioBucket = $XCTrioBuckets{$MomAlleles}{$DadAlleles}{$ChildAlleles};
                    } else {
                        $TrioBucket = $TrioBuckets{$MomAlleles}{$DadAlleles}{$ChildAlleles};
                    }
                    if (!defined($TrioBucket)) {
                        print "Couldn't find a bucket for pedigree $Ped, individual $Ind marker "
                          . $Loci[ $i + 1 ]
                          . ", [M]/[D]/[C] [$MomAlleles]/[$DadAlleles]/[$ChildAlleles], probably a Mendelian error!\n";
                        exit;
                    }

#			print "Pedigree $Ped / Marker ".$Loci[$i+1]." child $Ind (".$MomAlleles."-".$DadAlleles."-".$ChildAlleles.") gets bucket $TrioBucket\n";
                    # Add a child affection prefix and maybe a gender for XC analysis
                    if (defined($Directives{XC}) || $XC) {
                        push @bucketList, $TrioBucket . "-" . $Pedigrees{$Ped}{$Ind}{Sex} . $Pedigrees{$Ped}{$Ind}{Aff};
                    } else {
                        push @bucketList, $TrioBucket . "-" . $Pedigrees{$Ped}{$Ind}{Aff};
                    }
                }
            }
            my $PedBucket = $PAP . "/" . join("+", sort (@bucketList));
            $Buckets{ $Loci[ $i + 1 ] . "_" . $PedBucket }++;
            if (!defined($Templates{$PedBucket})) {
                $Templates{$PedBucket}{Ped}    = $Ped;
                $Templates{$PedBucket}{PedSeq} = sprintf("P%04d", $PedSeq++);
                $Templates{$PedBucket}{PairID} = $i;
            }
        }
    }

    if (scalar(keys %Buckets) == 0) {
        print "Bucketizing cannot reduce your pedigree count.\n";
    } else {
        print sprintf(
            "Bucketizing can reduce your evaluation count from %d to %d, or by %2d%%\n",
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
        open OUT, ">".$Prefix."_pedigrees.Dat";
    } else {
        open OUT, ">".$Prefix."_pedigrees.Pre";
    }

    # First the intact pedigrees
    for my $Ped (@Skippies) {
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            print OUT join(" ", ($Ped, $Ind, $Pedigrees{$Ped}{$Ind}{Dad}, $Pedigrees{$Ped}{$Ind}{Mom})) . " ";
            print OUT
              join(" ", ($Pedigrees{$Ped}{$Ind}{Kid1}, $Pedigrees{$Ped}{$Ind}{nPs}, $Pedigrees{$Ped}{$Ind}{nMs})) . " "
              if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
# &&& TBS CHANGE THIS TO LOOP OVER Loci &&&
            print OUT $Pedigrees{$Ped}{$Ind}{Aff} . " " . join(" ", @{ $Pedigrees{$Ped}{$Ind}{Mks} }) . "\n";
        }
    }

    # Next the template pedigrees
    print "Writing $Type pedigree\n";
    for my $PB (sort numericIsh keys %Templates) {
        my $Ped    = $Templates{$PB}{Ped};
        my $PairID = $Templates{$PB}{PairID};
        my $PedSeq = $Templates{$PB}{PedSeq};
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            print OUT join(" ", ($PedSeq, $Ind, $Pedigrees{$Ped}{$Ind}{Dad}, $Pedigrees{$Ped}{$Ind}{Mom})) . " ";
            print OUT
              join(" ", ($Pedigrees{$Ped}{$Ind}{Kid1}, $Pedigrees{$Ped}{$Ind}{nPs}, $Pedigrees{$Ped}{$Ind}{nMs})) . " "
              if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
            my $Pair = " " . $Pedigrees{$Ped}{$Ind}{Mks}[$PairID];
# &&& TBS CHANGE THIS TO USE REDUCED COUNT &&&
            print OUT $Pedigrees{$Ped}{$Ind}{Aff} . " " . join(" ", $Pair x $PairCount) . " ";
            print OUT "Ped: $PedSeq Per: $Ind";
            print OUT " # Template $PB\n";
        }
    }
    close OUT;

    if ($Type ne "POST") {
        system("makeped ".$Prefix."_pedigrees.Pre ".$Prefix."_pedigrees.Dat N");
    }

    # Finally the counts.
    open OUT, ">".$Prefix."_counts.Dat";

    print OUT "MARKER\t";
    for my $PB (sort numericIsh keys %Templates) {
        print OUT $Templates{$PB}{PedSeq} . "\t";
    }
    print OUT "\n";
    for my $i (0 .. $PairCount - 1) {
	next if (!$LociAttributes{$Loci[ $i + 1 ]}{Included});
        print OUT $Loci[ $i + 1 ] . "\t";
        for my $PB (sort numericIsh keys %Templates) {
            my $FB = $Loci[ $i + 1 ] . "_" . $PB;
            if (!defined($Buckets{$FB})) {
                print OUT "0\t";
            } else {
                print OUT $Buckets{$FB} . "\t";
            }
        }
        print OUT "\n";
    }
    close OUT;

    # If there was no configuration, create all of the supporting files
    if ($config) {
	print "Remember to modify your configuration file to specify the new pedigree and count files.\n";
	return;
    }

    open OUT, ">".$Prefix."_config.Dat";
    print OUT "PD ".$Prefix."_pedigrees.Dat\n";
    print OUT "DF ".$Prefix."_data.Dat\n";
    print OUT "MK ".$Prefix."_markers.Dat\n";
    print OUT "MP ".$Prefix."_map.Dat\n";
    print OUT "CC ".$Prefix."_counts.Dat\n";
    print OUT "HE ".$Prefix."_br.Out\n";
    print OUT "PF ".$Prefix."_ppl.Out\n";

    print OUT <<EOF;
PE
TP # Two-point analysis
Th 0 0.5 0.01
LD -1 1 0.1

# The rest is the standard analysis grid...
GF 0.001;0.01;0.1;0.3;0.5;0.8
DD 0.0 0.9 0.1
Dd 0.0 0.9 0.1
dd 0.0 0.9 0.1
DD 0.999
Dd 0.999
dd 0.999
DD >= Dd
Dd >= dd
DD != Dd; Dd != dd
AL 0.05 1 0.05

EOF
    print OUT "XC\n" if (defined($Directives{XC}) || $XC);
    close OUT;

    open OUT, ">".$Prefix."_data.Dat";
    for my $Name (@Loci) {
	next if (!$LociAttributes{$Name}{Included});
        print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
    }
    close OUT;

    open OUT, ">".$Prefix."_markers.Dat";
    for my $Name (@Loci) {
	next if (!$LociAttributes{$Name}{Included});
        if ($LociAttributes{$Name}{Type} eq "M") {
            print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
            print OUT "F";
            for my $Freq (@{ $LociAttributes{$Name}{Frequencies} }) {
                print OUT sprintf(" %.4f", $Freq);
            }
            print OUT "\n";
        }
    }
    close OUT;

    #
    open OUT, ">".$Prefix."_map.Dat";
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

    my $Type = shift();            # Pedigree type for writing
    my $Prefix = shift(); # Uniqifying (what a word!) prefix for files

    # Now write-out at least the pedigree and counts

    if ($Type eq "POST") {
        open OUT, ">".$Prefix."_pedigrees.Dat";
    } else {
        open OUT, ">".$Prefix."_pedigrees.Pre";
    }

    for my $Ped (sort numericIsh keys %Pedigrees ) {
        for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
            print OUT sprintf("%4s %3s %3s %3s ", $Ped, $Ind, $Pedigrees{$Ped}{$Ind}{Dad}, $Pedigrees{$Ped}{$Ind}{Mom});
            print OUT sprintf("%3s %3s %3s ", $Pedigrees{$Ped}{$Ind}{Kid1}, 
			      $Pedigrees{$Ped}{$Ind}{nPs}, $Pedigrees{$Ped}{$Ind}{nMs})
		if ($Type eq "POST");
            print OUT $Pedigrees{$Ped}{$Ind}{Sex} . " ";
            print OUT $Pedigrees{$Ped}{$Ind}{Prb} . " " if ($Type eq "POST");
	    print OUT $Pedigrees{$Ped}{$Ind}{Aff} . "  ";
            my @Pairs = @{ $Pedigrees{$Ped}{$Ind}{Mks} };
	    for my $i (0 .. $PairCount - 1) {
		next if (!$LociAttributes{$Loci[$i + 1]}{Included});
		next if ($LociAttributes{$Loci[$i + 1]}{Type} eq "T");
		next if ($LociAttributes{$Loci[$i + 1]}{Type} eq "A");
		print OUT $Pairs[$i] . "  "
	    }
            print OUT "\n";
        }
    }
    close OUT;

    if ($Type ne "POST") {
        system("makeped ".$Prefix."_pedigrees.Pre ".$Prefix."_pedigrees.Dat N");
    }

    open OUT, ">".$Prefix."_data.Dat";
    for my $Name (@Loci) {
	next if (!$LociAttributes{$Name}{Included});
        print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
    }
    close OUT;

    open OUT, ">".$Prefix."_markers.Dat";
    for my $Name (@Loci) {
	next if (!$LociAttributes{$Name}{Included});
        if ($LociAttributes{$Name}{Type} eq "M") {
            print OUT $LociAttributes{$Name}{Type} . " " . $Name . "\n";
            print OUT "F";
            for my $Freq (@{ $LociAttributes{$Name}{Frequencies} }) {
                print OUT sprintf(" %.4f", $Freq);
            }
            print OUT "\n";
        }
    }
    close OUT;

    #
    open OUT, ">".$Prefix."_map.Dat";
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
	$LociAttributes{$Name}{Included} = 1;
    }
    if (scalar(@include)) {
	# If we're doing inclusion, then turn them all off
	for my $Name (@Loci) {
	    $LociAttributes{$Name}{Included} = 0;
	}
	# Then turn on what is requested
	for my $Name (@include) {
	    $Name = sprintf("M%04d", $Name) if (!$config);
	    if (defined($LociAttributes{$Name}{Included})) {
		$LociAttributes{$Name}{Included} = 1;
	    } else {
		print "Marker \"$Name\" specified in -include list not found!\n";
	    }
	}
    }
    if (scalar(@exclude)) {
	# Turn off what is requested
	for my $Name (@exclude) {
	    $Name = sprintf("M%04d", $Name) if (!$config);
	    if (defined($LociAttributes{$Name}{Included})) {
		$LociAttributes{$Name}{Included} = 0;
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
sub kelvinConstraints {
    # General limitations
    if (defined($Directives{DK})) {
	# Integration (DK) analysis limitations
    } else {
	# Iterative (classic Kelvin) analysis limitations
    }
}

#####################################
# Verify command line parameters, check flags and do what the user asks.
#
my $Usage = <<EOF;

Usage "perl $0 [<flags>...] <input file>"

where <flags> are any of:

-config		The input file specified is a KELVIN configuration file. Otherwise
		it is assumed to be a pedigree file.
-pre		Pedigrees are in pre-MAKEPED format.
-post		Pedigrees are in post-MAKEPED format.
-noparents	Pedigrees are in pre-MAKEPED format with no columns for parents.
-XC		This is a sex-linked (X-chromosome) analysis
-bare		The pedigree file has only individual, affection status and marker
		allele pairs columns.
-nokelvin	Skip verification that kelvin can handle the analysis.
-loops		Check for consanguinity and marriage loops and print them if found.
-stats		Print statistics on the make-up of the pedigree(s).
-counts		Count genotypically identical pedigrees and print statistics.
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
PC_counts.Dat - counts of genotypically identical pedigrees.
PC_pedigrees.Pre or .Dat - pedigrees corresponding to counts, and any uncounted
        "pass-thru" pedigrees. If in pre-makeped format and makeped is found,
        it will be run to generate PC_pedigrees.Dat.

Output files generated by -write without the -config flag (<prefix> will be used
instead of "PC" if specified):
PC_data.Dat - pedigree companion data file indicating loci column usage.
PC_markers.Dat - marker data with allele freqs from unaffected (control) individuals.
PC_map.Dat - minimal dummy map file.
PC_config.Dat - template kelvin config file.

EOF

GetOptions(
	   'config'    => \$config,
	   'pre'       => \$pre,
	   'post'      => \$post,
	   'noparents' => \$noparents,
	   'XC'        => \$XC,
	   'bare'      => \$bare,
	   'nokelvin'  => \$nokelvin,
	   'loops'     => \$loops,
	   'stats'     => \$stats,
	   'count'     => \$count,
	   'include=s' => \@include,
	   'exclude=s' => \@exclude,
	   'split=i'   => \$split,
	   'write:s'   => \$write,
	   ) or die "Invalid command line parameters.";
if ($write ne "unspecified") {
    $WritePrefix = $write if ($write ne "");
    $write = 1;
} else {
    $write = 0;
}
@include = expand(/,/,join(',',@include)) if (@include);
@exclude = expand(/,/,join(',',@exclude)) if (@exclude);

die "Invalid number of arguments supplied.\n$Usage" if ($#ARGV < 0);
print "-config flag seen\n"                         if ($config);
print "-pre flag seen\n"                            if ($pre);
print "-post flag seen\n"                           if ($post);
print "-noparents flag seen\n"                      if ($noparents);
print "-XC flag seen\n"                             if ($XC);
print "-bare flag seen\n"                           if ($bare);
print "-nokelvin flag seen\n"                       if ($nokelvin);
print "-loops flag seen\n"                          if ($loops);
print "-stats flag seen\n"                          if ($stats);
print "-count flag seen\n"                          if ($count);
print "-include list of ".Dumper(\@include)." seen\n" if (@include);
print "-exclude list of ".Dumper(\@exclude)." seen\n" if (@exclude);
print "-split of $split seen\n"                     if ($split);
print "-write seen, using \"$WritePrefix\" prefix\n"   if ($write);
die "-pre -post and -bare are mutually exclusive flags."
  if ($pre + $post + $bare > 1);

if ($config) {
    my $ConfFile = shift;
    loadConf($ConfFile);
    if (defined($Directives{DF}[0])) { $companionFile = $Directives{DF}[0]; }
    loadCompanion($companionFile);
    if (defined($Directives{MP}[0])) { $mapFile = $Directives{MP}[0]; }
    loadMap($mapFile);
    if (defined($Directives{MK}[0])) { $markersFile = $Directives{MK}[0]; }
    loadMarkers($markersFile);
    if (defined($Directives{PD}[0])) { $pedFile = $Directives{PD}[0]; }
} else {
    $pedFile = shift;
}

my $pedFileType = assessPedigree($pedFile);
loadPedigree($pedFile, $pedFileType);

#print Dumper(\%Pedigrees);

if (!scalar(@Loci)) {

    # Make-up @Loci and %LociAttributes if we have to...
    deriveLociAndAttributes();
}

#print Dumper(\@Loci);
#print Dumper(\%LociAttributes);

checkRelations($pedFileType);
checkIntegrity();

if ($stats) {
    perfStats();
}

if ($loops) {
    if (!consanguinityLoop()) {
        marriageLoop();
    }
}

if ($nokelvin) {
    kelvinConstraints();
}

doMarkerInclusion();

if (defined($Directives{SA}) || defined($Directives{SS})) {
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
    my $IncludedMarkers = 0; my $SplitSet = 0;
    for my $i (1..$PairCount) {
	if (($LociAttributes{$Loci[$i]}{Included}) && (++$IncludedMarkers >= $split)) {
	    # Hit our limit, exclude all the rest
	    $SplitSet++;
	    for my $j (($i + 1)..$PairCount) {
		$LociAttributes{$Loci[$j]}{Included} = 0;
	    }
	    # Do the work
	    if ($count) {
		bucketizePedigrees($pedFileType,  $WritePrefix . $SplitSet);
	    } elsif ($write) {
		writeExpanded($pedFileType, $WritePrefix . $SplitSet);
	    }
	    # Redo the inclusion
	    doMarkerInclusion();
	    $IncludedMarkers = 0;
	}
    }
    if ($IncludedMarkers != 0) {
	# Do the rest
	if ($count) {
	    bucketizePedigrees($pedFileType,  $WritePrefix . $SplitSet);
	} elsif ($write) {
	    writeExpanded($pedFileType, $WritePrefix . $SplitSet);
	}
    }
}
exit;
