#!/usr/bin/perl -ws
# Strict doesn't co-exist with the -s option, so just uncomment as a diagnostic
#use strict;
use List::Util qw(sum);
use Data::Dumper;

# Usual header comments go here...

$| = 1;    # Force flush of output as printed.

# Sanctioned globals

# Permanent defaults
use constant AttributeMissing => 0;    # For marker alleles and Sex

# Defaults to be overridden by configuration file directives
my $pedFile          = "pedpost.dat";
my $companionFile    = "datafile.dat";
my $markersFile      = "markers.dat";    # Default input files
my $UnknownAffection = 0;
my $Unaffected       = 1;
my $Affected         = 2;                # Default affection indicators
my $UnknownPerson    = 0;                # Default unknown person indicator

# Read/calculated results
my %Pedigrees;                           # Pedigrees as loaded
my %Directives;                          # Directives as loaded
my $PairCount = 0;                       # Last pedigree count of marker pairs
my @Loci;                                # Ordered loci list from companion file
my %LociAttributes;                      # Loci attributes from companion and marker files

# Nuisances to fix
my @Depths;                              # Referenced recursively and I'm fuddled
my @Ancestors;                           # ditto
my $ShortestLoop = "";                   # Just batted around too much to monkey with right now.

my %KnownDirectives = (
    AL => \&NoAction,
    AS => \&dirAS,
    CT => \&dirCT,
    DA => \&NoAction,
    DD => \&NoAction,
    DF => \&NoAction,
    DT => \&NoAction,
    Dd => \&NoAction,
    GF => \&NoAction,
    HE => \&NoAction,
    MK => \&NoAction,
    MP => \&NoAction,
    PD => \&NoAction,
    PE => \&NoAction,
    QT => \&dirQT,
    SA => \&NoAction,
    SS => \&NoAction,
    TP => \&NoAction,
    Th => \&NoAction,
    UP => \&dirUP,
    dD => \&NoAction,
    dd => \&NoAction
);

#####################################
sub dirAS {
    $UnknownAffection = $Directives{AS}[0];
    $Unaffected       = $Directives{AS}[1];
    $Affected         = $Directives{AS}[2];
}

#####################################
sub dirCT {

    # If it's not already explicitly specified, set QT default for AS
    if (!defined($Directives{AS})) {
        $UnknownAffection = -99.99;
        $Unaffected       = -88.88;
        $Affected         = 88.88;
    }
}

#####################################
sub dirQT {

    # If it's not already explicitly specified, set QT default for AS
    if (!defined($Directives{AS})) {
        $UnknownAffection = -99.99;
        $Unaffected       = -88.88;
        $Affected         = 88.88;
    }
}

#####################################
sub dirUP {
    $UnknownPerson = $Directives{UP}[0];
}

#####################################
sub NoAction {
}

#####################################
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
    if ($Start == $End) {

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
                    if (   ($Pedigrees{$Pedily}{$Child}{Dad} == $Ind)
                        || ($Pedigrees{$Pedily}{$Child}{Mom} == $Ind)) {
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
    my ($Pedily, $Individual, $Depth);
    $Pedily     = shift();
    $Individual = shift();
    $Depth      = shift();

    $Depth++;
    my ($Dad, $Mom);
    $Dad = $Pedigrees{$Pedily}{$Individual}{Dad};
    $Mom = $Pedigrees{$Pedily}{$Individual}{Mom};
    if ($Dad != $UnknownPerson) {
        push @Depths,    $Depth;
        push @Ancestors, $Dad;
        listAncestors($Pedily, $Dad, $Depth);
    }
    if ($Mom != $UnknownPerson) {
        push @Depths,    $Depth;
        push @Ancestors, $Mom;
        listAncestors($Pedily, $Mom, $Depth);
    }
    return;
}

#####################################
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
sub numerically { $a <=> $b }

#####################################
sub consanguinityLoop() {

# One type of loop is two parents of a child, where the parents share common
# ancestry. Loop span is sum of two distances to common ancestor, so brother
# and sister parents would have a span of 2, because they both count back 1 to
# their mother (or father).
    my @Pairings = ();
    for my $Ped (keys %Pedigrees) {
        my %Seen = ();
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom != $UnknownPerson) && ($Dad != $UnknownPerson)) {
                push @Pairings, $Mom . "+" . $Dad unless $Seen{ $Mom . "+" . $Dad }++;
            }
        }
        for my $Pair (@Pairings) {
            my ($Mom, $Dad) = split /\+/, $Pair;
            my @Depths    = ();
            my @Ancestors = ();    # This could be better, but I'm cowardly
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
                    }
                }
            }
        }
    }
}

#####################################
sub marriageLoop() {

# The more extensive loop that cares about more than common ancestry is any path
# over or up from one parent back down or over to the other without retracing steps.
    for my $Ped (keys %Pedigrees) {
        my %Seen     = ();
        my @Pairings = ();
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom != $UnknownPerson) && ($Dad != $UnknownPerson)) {
                push @Pairings, $Mom . "+" . $Dad unless $Seen{ $Mom . "+" . $Dad }++;
            }
        }
        for my $Pair (@Pairings) {
            my ($Mom, $Dad) = split /\+/, $Pair;

            # Follow paths from Mom's Dad looking for Dad, having seen Mom
            my $DadsGrandPa = $Pedigrees{$Ped}{$Dad}{Dad};
            my $DadsGrandMa = $Pedigrees{$Ped}{$Dad}{Mom};
            if (   ($DadsGrandPa != $UnknownPerson)
                || ($DadsGrandMa != $UnknownPerson)) {

                my $MomsGrandPa = $Pedigrees{$Ped}{$Mom}{Dad};
                my $MomsGrandMa = $Pedigrees{$Ped}{$Mom}{Mom};
                if ($MomsGrandPa != $UnknownPerson) {
                    followPaths($Ped, $MomsGrandPa, $Dad, $Mom);
                }

                # Follow paths from Mom's Mom looking for Dad, having seen Mom
                if ($MomsGrandMa != $UnknownPerson) {
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
# Open and assess a pedigree file to determine processing required

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
# Open and read a post-makeped pedigree file producing a two-dimensional hash with
# the first being the family (pedigree) ID, and the second being the individual ID.
#
# Pedigrees{$Ped}{$Ind}
#
# ...at the moment, this is global and dynamic.

sub loadPedigree {
    my $File = shift();
    my $Type = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    %Pedigrees = ();
    my $MkC = 0;
    my $GtC;
    my ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Markers);

    while (<IN>) {
        $LineNo++;
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        if ($Type eq "POST") {
            ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Markers) = split /\s+/;
            die "Pedigree line $LineNo: proband must be numeric, not \"$Prb\"." if (!($Prb =~ /[0-9]/));
        } elsif ($Type eq "PRE") {
            ($Ped, $Ind, $Dad, $Mom, $Sex, $Aff, @Markers) = split /\s+/;
        } elsif ($Type eq "BARE") {
            ($Ind, $Aff, @Markers) = split /\s+/;
            $Ped = $Ind;
        } else {
            die "Unhandled pedigree type \"$Type\".";
        }

        # Validate everything we've got so far
        die "Pedigree line $LineNo: sex must be 1 or 2, not \"$Sex\"." if ((!defined($bare)) && !($Sex =~ /[12]/));
        die
          "Pedigree line $LineNo: affection status must be $UnknownAffection, $Unaffected, or $Affected, not \"$Aff\"."
          if (($Aff != $UnknownAffection) && ($Aff != $Unaffected) && ($Aff != $Affected));
        my ($OldFam, $OldInd, $Left);
        $MkC = 0;
        $GtC = 0;
        my @Pairs = ();
        for my $Marker (@Markers) {

            if ($Marker eq "Ped:") {
                $OldFam = $Markers[ $MkC + 1 ];
                $OldInd = $Markers[ $MkC + 3 ];
                last;
            }
            $MkC++;
            $GtC++ if ($Marker != AttributeMissing);
            if ($MkC & 1) {
                $Left = $Marker;
            } else {
                push @Pairs, $Left . " " . $Marker;
            }
        }
        if (!defined($bare)) {
            if (!defined($noparents)) {
                $Pedigrees{$Ped}{$Ind}{Dad} = $Dad;
                $Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
            }
            if (!defined($pre)) {
                $Pedigrees{$Ped}{$Ind}{Kid1} = $Kid1;
                $Pedigrees{$Ped}{$Ind}{nPs}  = $nPs;    # Next paternal sibling
                $Pedigrees{$Ped}{$Ind}{nMs}  = $nMs;    # Next maternal sibling
                $Pedigrees{$Ped}{$Ind}{Prb}  = $Prb;    # Proband
            }
            $Pedigrees{$Ped}{$Ind}{Sex} = $Sex;
        }
        $Pedigrees{$Ped}{$Ind}{Aff}    = $Aff;
        $Pedigrees{$Ped}{$Ind}{MkC}    = $MkC;          # Marker count (should always be the same)
        $Pedigrees{$Ped}{$Ind}{GtC}    = $GtC;          # Genotype count (how complete is genotyping)
        $Pedigrees{$Ped}{$Ind}{Mks}    = [@Pairs];
        $Pedigrees{$Ped}{$Ind}{OldFam} = $OldFam;
        $Pedigrees{$Ped}{$Ind}{OldInd} = $OldInd;
    }
    $PairCount = $MkC / 2;
    close IN;

    # Adopt single case/control individuals
    for my $Ped (sort keys %Pedigrees) {
	my $memberCount = scalar(keys %{ $Pedigrees{$Ped} });
	if ($memberCount == 1) {
	    for my $Ind (keys %{ $Pedigrees{$Ped} }) { # Yes, there's only one, but what individual ID?
		$Pedigrees{$Ped}{$Ind}{Prb} = 1;
		my $Dad = $Ind . "D";
		$Pedigrees{$Ped}{$Ind}{Dad} = $Dad;
		$Pedigrees{$Ped}{$Dad}{Sex} = 1;
		$Pedigrees{$Ped}{$Dad}{Dad} = $UnknownPerson;
		$Pedigrees{$Ped}{$Dad}{Mom} = $UnknownPerson;
		$Pedigrees{$Ped}{$Dad}{Aff} = $UnknownAffection;
		$Pedigrees{$Ped}{$Dad}{Prb} = 0;
		$Pedigrees{$Ped}{$Dad}{MkC} = 
		$Pedigrees{$Ped}{$Dad}{GtC} = 0;
		$Pedigrees{$Ped}{$Dad}{Mks} = [("0 0") x ($Pedigrees{$Ped}{$Ind}{MkC} /2 )];
		my $Mom = $Ind . "M";
		$Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
		$Pedigrees{$Ped}{$Mom}{Sex} = 2;
		$Pedigrees{$Ped}{$Mom}{Dad} = $UnknownPerson;
		$Pedigrees{$Ped}{$Mom}{Mom} = $UnknownPerson;
		$Pedigrees{$Ped}{$Mom}{Aff} = $UnknownAffection;
		$Pedigrees{$Ped}{$Mom}{Prb} = 0;
		$Pedigrees{$Ped}{$Mom}{MkC} = 
		$Pedigrees{$Ped}{$Mom}{GtC} = 0;
		$Pedigrees{$Ped}{$Mom}{Mks} = [("0 0") x ($Pedigrees{$Ped}{$Ind}{MkC} / 2)];
	    }
	}
    }
}

#####################################
# Derive a list of loci and marker alleles and frequencies from the pedigree itself

sub deriveLociAndAttributes {
    @Loci                        = ();
    %LociAttributes              = ();
    $Loci[0]                     = "Trait";
    $LociAttributes{Trait}{Type} = "T";
    for my $i (0 .. $PairCount - 1) {
        $Loci[$i] = sprintf("M%03d", $i);
        $LociAttributes{ $Loci[$i] }{Type} = "M";
        my @HaploCounts = (0); # List works because they must be numeric (relative position of frequency in marker file)
        for my $Ped (keys %Pedigrees) {
            for my $Ind (keys %{ $Pedigrees{$Ped} }) {
                my ($Left, $Right) = split /\s+/, $Pedigrees{$Ped}{$Ind}{Mks}[$i];
                $HaploCounts[$Left]++  if ($Left != AttributeMissing);
                $HaploCounts[$Right]++ if ($Right != AttributeMissing);
            }
        }

#	print "Looking at pair $i and seeing ".Dumper(\@HaploCounts)."\n";
        my $PopSize = sum @HaploCounts;
        my @Tokens  = ();
        for $i (1 .. scalar(@HaploCounts) - 1) {
            push @Tokens, $HaploCounts[$i] / $PopSize;
        }
        $LociAttributes{ $Loci[$i] }{Frequencies} = [@Tokens];
    }
}

#####################################
# Open and read the marker description companion file producing an ordered list of
# loci names and a hash for named loci attributes.

sub loadCompanion {
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    @Loci           = ();
    %LociAttributes = ();
    while (<IN>) {
        $LineNo++;
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        s/^\s*//g;       # Trim leading whitespace
        my ($Type, $Name) = split /\s+/;
        die "Unknown locus type \"$Type\" at line $LineNo in marker description companion file $File\n"
          if (($Type ne "T") and ($Type ne "M"));
        push @Loci, $Name;
        $LociAttributes{$Name}{Type} = $Type;
    }
    close IN;
}

#####################################
# Open and read the marker file to flesh-out the %LociAttributes hash

sub loadMarkers {
    my $File = shift();
    die "$File is not a file." if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    my $Name   = "";
    while (<IN>) {
        $LineNo++;
        s/\s*\#.*//g;    # Trim comments
        next if (/^$/);  # Drop empty lines
        my @Tokens     = split /\s+/;
        my $RecordType = shift @Tokens;
        if ($RecordType eq "M") {
            $Name = shift @Tokens;
        } elsif ($RecordType eq "F") {
            $LociAttributes{$Name}{Frequencies} = [@Tokens];
        } else {
            die "Unknown record type \"$RecordType\" at line $LineNo in marker file $File\n";
        }
    }
    close IN;
}

#####################################
# Check relations and genders (parent, siblings present, not self-parenting, correct genders)

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

# Start discovering statistics that might affect performance.
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
            if ($Pedigrees{$Ped}{$Ind}{Sex} == AttributeMissing)  { $PedCat{$Ped}{MsgSex}++; }
        }
        push @PedSiz, scalar(keys(%{ $Pedigrees{$Ped} }));
        push @PedNuk, scalar(keys(%Seen));
    }

    @PedSiz = @PedSiz;
    print "Number of pedigrees:" . @PedSiz . "\n";
    $Avg = sprintf("%.2f", (sum @PedSiz) / @PedSiz);
    @PedNuk = @PedNuk;
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
        @Counts = @Counts;
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
        @Ancestors = @Ancestors;
        $Avg = sprintf("%.2f", (sum @Ancestors) / @Ancestors);
        print "Pedigree $Ped Depth: max:" . $Ancestors[-1] . " med:" . $Ancestors[ $#Ancestors / 2 ] . " avg: $Avg\n";
    }
}

#####################################
#
# Adapted from a trios-only version by John Burian.
#
# Leverages pedigree inheritance patterns to reduce computational
# complexity. Multiple pedigrees that fit a single inheritance pattern
# are converted to a single pedigree and a weight for the likelihood poly.
# Currently works only for two-point analysis because it would be a real
# trick for there to be many pedigrees where multiple markers all fit the
# same inheritance pattern.
#
# Reads a map and pedigree file including affected sibling pair (ASP), trio
# (nuclear) and generates a new pedigree file and count file with redundant
# pedigrees counted by marker. Unique pedigrees are retained with a count of 1.
#
# Trios, i.e. 3-member families with an affected child and parents treated
# as if unaffected, fall into one of 30 categories.
#
# ASPs are expanded trios, i.e. 4-member pedigrees with two affected children
# and both parents in any affectation category, which will fall first into
# one of the 3-member trio categories for the first child, and then a second
# 3-member category with the same parentage for the second child, e.g. 30 and
# then 18. This is used to identify a unique category by ordering the 3-member
# numbers, e.g. 18-30, because a family falling into 30 and then 18 is the
# same as one falling into 18 and then 30. While this extends to any number
# of children of a single-generation family, we don't have such a thing as
# "quads" or "affected sibling triples" just yet. Generalizing this to deal with
# any single-generation family inheritance pattern would require incorporating
# affection status for all individuals. Do-able, but probably not worth the
# effort.
#
# Case/control is different because the parents are assumed to be unknown, and
# are generated if not present.
#
sub bucketizePedigrees {

    my %TrioBuckets = (
        '0 0' => {
            '0 0' => {
                '0 0' => 'T30',    # 30  0 0  0 0  0 0
                '1 1' => 'T18',    # 18  0 0  0 0  1 1
                '1 2' => 'T19',    # 19  0 0  0 0  1 2
                '2 2' => 'T20'     # 20  0 0  0 0  2 2
            },
            '1 1' => {
                '0 0' => 'T24',    # 24  0 0  1 1  0 0
                '1 1' => 'T25',    #  5  0 0  1 1  1 1
                '1 2' => 'T06'     #  6  0 0  1 1  1 2
            },
            '1 2' => {
                '0 0' => 'T27',    # 27  0 0  1 2  0 0
                '1 1' => 'T12',    # 12  0 0  1 2  1 1
                '1 2' => 'T13',    # 13  0 0  1 2  1 2
                '2 2' => 'T14'     # 14  0 0  1 2  2 2
            },
            '2 2' => {
                '0 0' => 'T29',    # 29  0 0  2 2  0 0
                '1 2' => 'T16',    # 16  0 0  2 2  1 2
                '2 2' => 'T17'     # 17  0 0  2 2  2 2
            }
        },
        '1 1' => {
            '0 0' => {
                '0 0' => 'T24',    # 24  1 1  0 0  0 0
                '1 1' => 'T05',    #  5  1 1  0 0  1 1
                '1 2' => 'T06'     #  6  1 1  0 0  1 2
            },
            '1 1' => {
                '0 0' => 'T21',    # 21  1 1  1 1  0 0
                '1 1' => 'T01'     #  1  1 1  1 1  1 1
            },
            '1 2' => {
                '0 0' => 'T22',    # 22  1 1  1 2  0 0
                '1 1' => 'T02',    #  2  1 1  1 2  1 1
                '1 2' => 'T03'     #  3  1 1  1 2  1 2
            },
            '2 2' => {
                '0 0' => 'T23',    # 23  1 1  2 2  0 0
                '1 2' => 'T04'     #  4  1 1  2 2  1 2
            }
        },
        '1 2' => {
            '0 0' => {
                '0 0' => 'T27',    # 27  1 2  0 0  0 0
                '1 1' => 'T12',    # 12  1 2  0 0  1 1
                '1 2' => 'T13',    # 13  1 2  0 0  1 2
                '2 2' => 'T14'     # 14  1 2  0 0  2 2
            },
            '1 1' => {
                '0 0' => 'T22',    # 22  1 2  1 1  0 0
                '1 1' => 'T02',    #  2  1 2  1 1  1 1
                '1 2' => 'T03'     #  3  1 2  1 1  1 2
            },
            '1 2' => {
                '0 0' => 'T25',    # 25  1 2  1 2  0 0
                '1 1' => 'T07',    #  7  1 2  1 2  1 1
                '1 2' => 'T08',    #  8  1 2  1 2  1 2
                '2 2' => 'T09'     #  9  1 2  1 2  2 2
            },
            '2 2' => {
                '0 0' => 'T26',    # 26  1 2  2 2  0 0
                '1 2' => 'T10',    # 10  1 2  2 2  1 2
                '2 2' => 'T11'     # 11  1 2  2 2  2 2
            }
        },
        '2 2' => {
            '0 0' => {
                '0 0' => 'T29',    # 29  2 2  0 0  0 0
                '1 2' => 'T16',    # 16  2 2  0 0  1 2
                '2 2' => 'T17'     # 17  2 2  0 0  2 2
            },
            '1 1' => {
                '0 0' => 'T23',    # 23  2 2  1 1  0 0
                '1 2' => 'T04'     #  4  2 2  1 1  1 2
            },
            '1 2' => {
                '0 0' => 'T26',    # 26  2 2  1 2  0 0
                '1 2' => 'T10',    # 10  2 2  1 2  1 2
                '2 2' => 'T11'     # 11  2 2  1 2  2 2
            },
            '2 2' => {
                '0 0' => 'T28',    # 28  2 2  2 2  0 0
                '2 2' => 'T15'     # 15  2 2  2 2  2 2
            }
        }
    );

    # Verify that this is a 2pt analysis (default, so look for multipoint directives)
    die "Generation of counts not permitted for a multipoint analysis.\n"
      if (defined($Directives{SA}) || defined($Directives{SS}));

    # Verify that all markers are present and only biallelic...
    for my $Marker (@Loci) {
        next if ($LociAttributes{$Marker}{Type} eq "T");
        die "No allele information found for marker $Marker for count generation.\n"
          if (!defined($LociAttributes{$Marker}{Frequencies}));
        die "Marker $Marker not biallelic, not permitted for count generation.\n"
          if (scalar(@{ $LociAttributes{$Marker}{Frequencies} }) != 2);
    }

    my %Buckets = ();

    # Look at each family...
    for my $Ped (sort keys %Pedigrees) {

	my $memberCount = scalar(keys %{ $Pedigrees{$Ped} });
        # Qualify the family for inclusion in trio buckets...
        if (($memberCount != 3) && ($memberCount != 4)) {
            print "Copy $memberCount-individual pedigree $Ped directly\n";
            next;
        }

	# Build a parental affection prefix so we can do more than expected.
	my $PAP = "";
        for my $Ind (sort keys %{ $Pedigrees{$Ped} }) {
            if ($Pedigrees{$Ped}{$Ind}{Dad} eq $UnknownPerson) {
		if ($Pedigrees{$Ped}{$Ind}{Sex} == 1) {
		    $PAP = $Pedigrees{$Ped}{$Ind}{Aff}.$PAP ;
		} else {
		    $PAP = $PAP.$Pedigrees{$Ped}{$Ind}{Aff};
		}
	    }
	}

        # Get the family genotype bucket for each marker pair
        for my $i (0 .. $PairCount - 1) {
            my @bucketList = ();

            # Get a trio bucket for each child in the family
            for my $Ind (sort keys %{ $Pedigrees{$Ped} }) {

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
                    my $TrioBucket = $TrioBuckets{$MomAlleles}{$DadAlleles}{$ChildAlleles};
                    if ($TrioBucket eq "") {
                        print "Couldn't find a bucket for pedigree $Ped, marker "
                          . $Loci[$i]
                          . ", M/D/C $MomAlleles/$DadAlleles/$ChildAlleles!\n";
                        exit;
                    }
#			print "Pedigree $Ped / Marker ".$Loci[$i]." child $Ind (".$MomAlleles."-".$DadAlleles."-".$ChildAlleles.") gets bucket $TrioBucket\n";
		    # Add a child affection prefix
                    push @bucketList, $TrioBucket . "+" . $Pedigrees{$Ped}{$Ind}{Aff};
                }
            }
            my $FullBucket = $PAP . "_" . join("-", sort (@bucketList)) . "_" . $Loci[$i];

#		print "Ped $Ped goes into full bucket $FullBucket for marker ".$Loci[$i]."\n";
            if (!defined($Buckets{$FullBucket})) {

#                    print "Defining bucket [$FullBucket]\n";
                $Buckets{$FullBucket} = 1;
            } else {

#		    print "Bumping bucket [$FullBucket]\n";
                $Buckets{$FullBucket}++;
            }
        }
    }
    for my $Bucket (sort keys(%Buckets)) {
        print "Bucket $Bucket has " . $Buckets{$Bucket} . " entries\n";
    }
}

# Verify command line parameters
my $Usage = <<EOF;

Usage "perl $0 [<flags>...] <input file>"

where <flags> are any of:

-config		The input file specified is a KELVIN configuration file. Otherwise
		it is assumed to be a pedigree file.
-pre		Pedigrees are in pre-MAKEPED format.
-post		Pedigrees are in post-MAKEPED format.
-noparents	Pedigrees are in pre-MAKEPED format with no columns for parents.
-bare		The pedigree file has only affection status and marker allele pairs
-counts		Generate new pedigree file and counts

The input file will be read and analyzed. If it is a configuration file,
the input files it references will be read and analyzed as well. 

Validates pedigrees and other configuration files, finds loops, expands
case-control individuals, generates count files for genotypically 
identical pedigrees.

If counts are requested, a new pedigree file (in the same PRE/POST format 
as the input) will be generated along with a count file. If only a pedigree
file was provided, dummy configuration, pedigree companion, marker and map 
files will be generated for a two-point LD run with marker allele frequencies 
computed from the controls?! If a configuration file was provided, the marker 
names from the companion and marker files it specifies will be used.

EOF

die "Invalid number of arguments supplied.\n$Usage" if ($#ARGV < 0);
print "-config flag seen\n"                         if ($config);
print "-pre flag seen\n"                            if ($pre);
print "-post flag seen\n"                           if ($post);
print "-noparents flag seen\n"                      if ($noparents);
print "-bare flag seen\n"                           if ($bare);
print "-counts flag seen\n"                         if ($counts);
die "-pre -post and -bare are mutually exclusive flags."
  if (defined($pre) + defined($post) + defined($bare) > 1);

if ($config) {
    my $ConfFile = shift;
    print "Processing configuration file $ConfFile\n";
    loadConf($ConfFile);
    if (defined($Directives{DF}[0])) { $companionFile = $Directives{DF}[0]; }
    loadCompanion($companionFile);
    if (defined($Directives{MK}[0])) { $markersFile = $Directives{MK}[0]; }
    loadMarkers($markersFile);
    if (defined($Directives{PD}[0])) { $pedFile = $Directives{PD}[0]; }
} else {
    $pedFile = shift;
    print "Processing pedigree file $pedFile\n";
}

my $pedFileType = assessPedigree($pedFile);
print "Processing $pedFileType format pedigree file\n";
loadPedigree($pedFile, $pedFileType);

print Dumper(\%Pedigrees);

if (!scalar(@Loci)) {

    # Make-up @Loci and %LociAttributes if we have to...
    deriveLociAndAttributes();
}

#print Dumper(\@Loci);
#print Dumper(\%LociAttributes);

checkRelations($pedFileType);
#    perfStats();
if (!defined($bare)) {
    consanguinityLoop();
    marriageLoop();
}

if ($counts) {
    bucketizePedigrees();
}
exit;
