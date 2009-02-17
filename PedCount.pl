#!/usr/bin/perl -w
use List::Util qw(sum);
use strict;

# Usual header comments go here...

$| = 1;    # Force flush of output as printed.

use constant FALSE => 0;
use constant TRUE => 1;

# Defaults to be overridden by configuration file directives
use constant AttributeMissing => 0;
use constant DefaultUnknownAffection => 0;
use constant DefaultUnaffected   => 1;
use constant DefaultAffected => 2;

# Sanctioned globals
my $Usage = "Usage \"perl $0 <configuration file, like kelvin.conf>\"\n";
my @Depths; # Referenced recursively and I'm fuddled
my @Ancestors; # ditto
my $ShortestLoop = "";
my ($UnknownAffection, $Unaffected, $Affected);
my %Pedigrees; # Pedigrees as loaded
my %Directives; # Directives as loaded
my $PairCount = 0; # Last pedigree count of marker pairs
my @Loci; # Ordered loci list from companion file
my %LociAttributes; # Loci attributes from companion and marker files
my %TrioBuckets = 
    ( ' 0 0' => { ' 0 0' => { ' 0 0' => 'T30',  # 30  0 0  0 0  0 0
			      ' 1 1' => 'T18',  # 18  0 0  0 0  1 1
			      ' 1 2' => 'T19',  # 19  0 0  0 0  1 2
			      ' 2 2' => 'T20'   # 20  0 0  0 0  2 2
			      },
		  ' 1 1' => { ' 0 0' => 'T24',  # 24  0 0  1 1  0 0
			      ' 1 1' => 'T25',  #  5  0 0  1 1  1 1
			      ' 1 2' => 'T06'   #  6  0 0  1 1  1 2
			      },
		  ' 1 2' => { ' 0 0' => 'T27',  # 27  0 0  1 2  0 0
			      ' 1 1' => 'T12',  # 12  0 0  1 2  1 1
			      ' 1 2' => 'T13',  # 13  0 0  1 2  1 2
			      ' 2 2' => 'T14'   # 14  0 0  1 2  2 2
		  },
		  ' 2 2' => { ' 0 0' => 'T29',  # 29  0 0  2 2  0 0
			      ' 1 2' => 'T16',  # 16  0 0  2 2  1 2
			      ' 2 2' => 'T17'   # 17  0 0  2 2  2 2
		  }
	      },
      ' 1 1' => { ' 0 0' => { ' 0 0' => 'T24',  # 24  1 1  0 0  0 0
                              ' 1 1' => 'T05',  #  5  1 1  0 0  1 1
                              ' 1 2' => 'T06'   #  6  1 1  0 0  1 2
			      },
		  ' 1 1' => { ' 0 0' => 'T21',  # 21  1 1  1 1  0 0
                              ' 1 1' => 'T01'   #  1  1 1  1 1  1 1
			       },
		  ' 1 2' => { ' 0 0' => 'T22',  # 22  1 1  1 2  0 0
                              ' 1 1' => 'T02',  #  2  1 1  1 2  1 1
                              ' 1 2' => 'T03'   #  3  1 1  1 2  1 2
			       },
		  ' 2 2' => { ' 0 0' => 'T23',  # 23  1 1  2 2  0 0
                              ' 1 2' => 'T04'   #  4  1 1  2 2  1 2
			       }
	      },
      ' 1 2' => { ' 0 0' => { ' 0 0' => 'T27',  # 27  1 2  0 0  0 0
                              ' 1 1' => 'T12',  # 12  1 2  0 0  1 1
                              ' 1 2' => 'T13',  # 13  1 2  0 0  1 2
                              ' 2 2' => 'T14'   # 14  1 2  0 0  2 2
			       },
		  ' 1 1' => { ' 0 0' => 'T22',  # 22  1 2  1 1  0 0
                              ' 1 1' => 'T02',  #  2  1 2  1 1  1 1
                              ' 1 2' => 'T03'   #  3  1 2  1 1  1 2
			       },
		  ' 1 2' => { ' 0 0' => 'T25',  # 25  1 2  1 2  0 0
                              ' 1 1' => 'T07',  #  7  1 2  1 2  1 1
                              ' 1 2' => 'T08',  #  8  1 2  1 2  1 2
                              ' 2 2' => 'T09'   #  9  1 2  1 2  2 2
			       },
		  ' 2 2' => { ' 0 0' => 'T26',  # 26  1 2  2 2  0 0
                              ' 1 2' => 'T10',  # 10  1 2  2 2  1 2
                              ' 2 2' => 'T11'   # 11  1 2  2 2  2 2
			       }
	      },
      ' 2 2' => { ' 0 0' => { ' 0 0' => 'T29',  # 29  2 2  0 0  0 0
                              ' 1 2' => 'T16',  # 16  2 2  0 0  1 2
                              ' 2 2' => 'T17'   # 17  2 2  0 0  2 2
			       },
		  ' 1 1' => { ' 0 0' => 'T23',  # 23  2 2  1 1  0 0
                              ' 1 2' => 'T04'   #  4  2 2  1 1  1 2
			       },
		  ' 1 2' => { ' 0 0' => 'T26',  # 26  2 2  1 2  0 0
                              ' 1 2' => 'T10',  # 10  2 2  1 2  1 2
                              ' 2 2' => 'T11'   # 11  2 2  1 2  2 2
			       },
		  ' 2 2' => { ' 0 0' => 'T28',  # 28  2 2  2 2  0 0
                              ' 2 2' => 'T15'   # 15  2 2  2 2  2 2
			  }
	      }
      );


#####################################
sub followPaths {
    my ($Pedily, $Start, $End, $Seen, $i, $Ind);
    $Pedily = shift();
    $Start  = shift();
    $End    = shift();
    $Seen   = shift();

    if ((++$i % 10000) == 0) {

#	print "$i: Fam $Pedily: $Seen -> $Start\n";
    }

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
    if ($Dad != AttributeMissing) {
        followPaths ($Pedily, $Dad, $End, $Seen);
    }
    $Mom = $Pedigrees{$Pedily}{$Start}{Mom};
    if ($Mom != AttributeMissing) {
        followPaths ($Pedily, $Mom, $End, $Seen);
    }

    # Try down (for every individual who has $Start as a parent)
    for $Ind (keys %{ $Pedigrees{$Pedily} }) {
        if (   ($Pedigrees{$Pedily}{$Ind}{Dad} == $Start)
            || ($Pedigrees{$Pedily}{$Ind}{Mom} == $Start)) {
            followPaths ($Pedily, $Ind, $End, $Seen);
        }
    }

    # Try my spouses (for every individual I share parenting with)
    my @Spouses    = ();
    my %SpouseSeen = ();
    for $Ind (keys %{ $Pedigrees{$Pedily} }) {
        if ($Start == $Pedigrees{$Pedily}{$Ind}{Mom}) {
            push @Spouses, $Pedigrees{$Pedily}{$Ind}{Dad} unless $SpouseSeen{ $Pedigrees{$Pedily}{$Ind}{Dad} }++;
        }
        if ($Start == $Pedigrees{$Pedily}{$Ind}{Dad}) {
            push @Spouses, $Pedigrees{$Pedily}{$Ind}{Mom} unless $SpouseSeen{ $Pedigrees{$Pedily}{$Ind}{Mom} }++;
        }
    }
    for my $Spouse (@Spouses) {
        followPaths ($Pedily, $Spouse, $End, $Seen);
    }

    # Try over (for every individual with a same parent as me but me)
    # but if both parents are the same and the sibling has no kids, don't bother.
    for $Ind (keys %{ $Pedigrees{$Pedily} }) {
        my $YourDad = $Pedigrees{$Pedily}{$Ind}{Dad};
        my $MyDad   = $Pedigrees{$Pedily}{$Start}{Dad};
        my $YourMom = $Pedigrees{$Pedily}{$Ind}{Mom};
        my $MyMom   = $Pedigrees{$Pedily}{$Start}{Mom};
        if (   (($YourMom != AttributeMissing) && ($YourMom == $MyMom))
            || (($YourDad != AttributeMissing) && ($YourDad == $MyDad))) {
            if (($YourMom == $MyMom) && ($YourDad == $MyDad)) {
                for my $Child (keys %{ $Pedigrees{$Pedily} }) {
                    if (   ($Pedigrees{$Pedily}{$Child}{Dad} == $Ind)
                        || ($Pedigrees{$Pedily}{$Child}{Mom} == $Ind)) {
                        followPaths ($Pedily, $Ind, $End, $Seen);
                        last;
                    }
                }
            } else {
                followPaths ($Pedily, $Ind, $End, $Seen);
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
    $Mom =  $Pedigrees{$Pedily}{$Individual}{Mom};
    if ($Dad != AttributeMissing) {
        push @Depths,    $Depth;
        push @Ancestors, $Dad;
        listAncestors ($Pedily, $Dad, $Depth);
    }
    if ($Mom != AttributeMissing) {
        push @Depths,    $Depth;
        push @Ancestors, $Mom;
        listAncestors ($Pedily, $Mom, $Depth);
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
    if ($Dad == AttributeMissing) {
        $DadCount = 0;
    } else {
        $DadCount = 1 + countAncestors ($Pedily, $Dad, $Count);
    }
    if ($Mom == AttributeMissing) {
        $MomCount = 0;
    } else {
        $MomCount = 1 + countAncestors ($Pedily, $Mom, $Count);
    }
    return ($DadCount + $MomCount + $Count);
}

#####################################
sub numerically { $a <=> $b }

#####################################
sub dump() {
# Dump the whole kit and kaboodle
    for my $Ped (sort keys %Pedigrees) {
        printf "Family $Ped\n";
        for my $Ind (sort keys %{ $Pedigrees{$Ped} }) {
            print "\t$Ind:";
            for my $Attribute (@{ $Pedigrees{$Ped}{$Ind} }) {
                printf(" $Attribute");
            }
            print "\n";
        }
    }
}


#####################################
sub sharedAncestry() {

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
            if (($Mom != AttributeMissing) && ($Dad != AttributeMissing)) {
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
            @Depths       = ();
            @Ancestors    = ();
            listAncestors($Ped, $Dad, 0);
            my @DadDepths    = @Depths;
            my @DadAncestors = @Ancestors;

            for my $i (0 .. $#MomDepths) {
                for my $j (0 .. $#DadDepths) {
                    if ($MomAncestors[$i] == $DadAncestors[$j]) {
                        my $LoopSize = $MomDepths[$i] + $DadDepths[$j];
                        print "Family $Ped loop of size $LoopSize at ancestor "
                          . $MomAncestors[$i]
                          . " for pair $Mom / $Dad\n";
                    }
                }
            }
        }
    }
}

#####################################
sub simpleLoop() {

# The more extensive loop that cares about more than common ancestry is any path
# over or up from one parent back down or over to the other without retracing steps.
    for my $Ped (keys %Pedigrees) {
        my %Seen     = ();
        my @Pairings = ();
        for my $Ind (keys %{ $Pedigrees{$Ped} }) {
            my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom != AttributeMissing) && ($Dad != AttributeMissing)) {
                push @Pairings, $Mom . "+" . $Dad unless $Seen{ $Mom . "+" . $Dad }++;
            }
        }
        for my $Pair (@Pairings) {
            my ($Mom, $Dad) = split /\+/, $Pair;

            # Follow paths from Mom's Dad looking for Dad, having seen Mom
            my $DadsGrandPa = $Pedigrees{$Ped}{$Dad}{Dad};
            my $DadsGrandMa = $Pedigrees{$Ped}{$Dad}{Mom};
            if (   ($DadsGrandPa != AttributeMissing)
                || ($DadsGrandMa != AttributeMissing)) {

                my $MomsGrandPa = $Pedigrees{$Ped}{$Mom}{Dad};
                my $MomsGrandMa = $Pedigrees{$Ped}{$Mom}{Mom};
                if ($MomsGrandPa != AttributeMissing) {
                    followPaths($Ped, $MomsGrandPa, $Dad, $Mom);
                }

                # Follow paths from Mom's Mom looking for Dad, having seen Mom
                if ($MomsGrandMa != AttributeMissing) {
                    followPaths($Ped, $MomsGrandMa, $Dad, $Mom);
                }
                if ($ShortestLoop ne "") {
                    print "Family $Ped loop for pair $Pair as [$ShortestLoop]\n";
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
    die "$File is not a file. $Usage" if (!-f $File);
    open IN, "<$File" || die "Cannot open $File\n";
    %Directives = ();
    while (<IN>) {
	s/\s*\#.*//g; # Trim comments
	next if (/^$/); # Drop empty lines
	my @Parameters = split;
	my $Directive = shift @Parameters;
	$Directives{$Directive} = \@Parameters;
    }
    if (defined($Directives{AS})) {
	$UnknownAffection = $Directives{AS}[0];
	$Unaffected = $Directives{AS}[1];
	$Affected = $Directives{AS}[2];
    } else {
	$UnknownAffection = DefaultUnknownAffection;
	$Unaffected = DefaultUnaffected;
	$Affected = DefaultAffected;
    }
}

#####################################
# Open and assess a pedigree file to determine processing required

sub assessPedigree {
    my $File = shift();
    die "$File is not a file. $Usage" if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    while (<IN>) {
	s/\s*\#.*//g; # Trim comments
	next if (/^$/); # Drop empty lines
	my @Columns = split /\s+/;
	# See if the number of columns makes any kind of sense...
	# Subtract 2*Loci - 1 (for single trait column) and 5 common columns
	my $Slack = scalar(@Columns) - ((scalar(@Loci) * 2) - 1) - 5;
	if ($Slack == 8) {
	    print "POST\n";
	} elsif ($Slack == 4) {
	    print "POST with old tail cut off\n";
	} elsif ($Slack == 0) {
	    print "PRE\n";
	} elsif ($Slack == -2) {
	    print "No-parent PRE (CC)\n";
	}
	close IN;
	return;
    }
}

#####################################
# Open and read a post-makeped pedigree file producing a two-dimensional hash with
# the first being the family (pedigree) ID, and the second being the individual ID.
#
# Pedigrees{$Ped}{$Ind}
#
# ...at the moment, this is global and dynamic.

sub loadPostPedigree {
    my $File = shift();
    die "$File is not a file. $Usage" if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    %Pedigrees = ();
    my $MkC;
    while (<IN>) {
	s/\s*\#.*//g; # Trim comments
	next if (/^$/); # Drop empty lines
	my ($Ped, $Ind, $Dad, $Mom, $Kid1, $nPs, $nMs, $Sex, $Prb, $Aff, @Markers) = split /\s+/;
	my ($OldFam, $OldInd, $Left);
	$MkC = 0;
	my @Pairs = ();
	for my $Marker (@Markers) {
	    if ($Marker eq "Ped:") {
		$OldFam = $Markers[$MkC+1];
		$OldInd = $Markers[$MkC+3];
		last;
	    }
	    $MkC++;
	    if ($MkC & 1) {
		$Left = $Marker;
	    } else {
		push @Pairs, " ".$Left." ".$Marker;
	    }
	}
	$Pedigrees{$Ped}{$Ind}{Dad} = $Dad;
	$Pedigrees{$Ped}{$Ind}{Mom} = $Mom;
	$Pedigrees{$Ped}{$Ind}{Kid1} = $Kid1;
	$Pedigrees{$Ped}{$Ind}{nPs} = $nPs;
	$Pedigrees{$Ped}{$Ind}{nMs} = $nMs;
	$Pedigrees{$Ped}{$Ind}{Sex} = $Sex;
	$Pedigrees{$Ped}{$Ind}{Prb} = $Prb;
	$Pedigrees{$Ped}{$Ind}{Aff} = $Aff;
	$Pedigrees{$Ped}{$Ind}{MkC} = $MkC;
	$Pedigrees{$Ped}{$Ind}{Mks} = [@Pairs];
	$Pedigrees{$Ped}{$Ind}{OldFam} = $OldFam;
	$Pedigrees{$Ped}{$Ind}{OldInd} = $OldInd;
    }
    $PairCount = ($MkC / 2) - 1;
    close IN;
}

#####################################
# Open and read the marker description companion file producing an ordered list of
# loci names and a hash for named loci attributes.

sub loadCompanion {
    my $File = shift();
    die "$File is not a file. $Usage" if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    @Loci = (); %LociAttributes = ();
    while (<IN>) {
	$LineNo++;
	s/\s*\#.*//g; # Trim comments
	next if (/^$/); # Drop empty lines
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
    die "$File is not a file. $Usage" if (!-f $File);
    open IN, "<$File" || die "Cannot open file $File\n";
    my $LineNo = 0;
    my $Name = "";
    while (<IN>) {
	$LineNo++;
	s/\s*\#.*//g; # Trim comments
	next if (/^$/); # Drop empty lines
	my @Tokens = split /\s+/;
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

    for my $Ped (sort keys %Pedigrees) {
	for my $Ind (keys %{ $Pedigrees{$Ped} }) {
	    my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            if ($Dad != AttributeMissing) {
		die "In family $Ped, $Ind is own father!\n"                        if ($Dad == $Ind);
		die "In family $Ped, Father $Dad not found for individual $Ind!\n" if (!defined($Pedigrees{$Ped}{$Dad}));
		die "In family $Ped, Father $Dad is not male!\n"                   if ($Pedigrees{$Ped}{$Dad}{Sex} != 1);
            }
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if ($Mom != AttributeMissing) {
		die "In family $Ped, $Ind is own mother!\n"                        if ($Mom == $Ind);
		die "In family $Ped, Mother $Mom not found for individual $Ind!\n" if (!defined($Pedigrees{$Ped}{$Mom}));
		die "In family $Ped, Mother $Mom is not female!\n"                 if ($Pedigrees{$Ped}{$Mom}{Sex} != 2);
            }
            my $Kid1 = $Pedigrees{$Ped}{$Ind}{Kid1};
            if ($Kid1 != AttributeMissing) {
		die "In family $Ped, $Ind is own sibling!\n" if ($Kid1 == $Ind);
                die "In family $Ped, first child $Kid1 missing for individual $Ind!\n" if (!defined($Pedigrees{$Ped}{$Kid1}));
            }
        }
    }
}

# Start discovering statistics that might affect performance.
sub perfStats {

    my %CountsLabels = (MsgMkr => "Missing Markers",
			MsgSib => "Missing Next Sibling(s)",
			MsgAff => "Missing Affected Status",
			MsgSex => "Missing Sex",
			Fdrs => "Founders",
			NoFdrs => "Non-Founders",
			FdrsMM => "Founders Without Markers", 
			MsgPPM => "Individuals Missing Parental Markers");

    my %PedMsg = ();
    my @PedSiz = ();
    my @PedNuk = ();
    my $Avg;
    for my $Ped (keys %Pedigrees) {
	push @PedSiz, keys(%{ $Pedigrees{$Ped} }) + 0;
	my %Seen = ();
	for my $Ind (sort keys %{ $Pedigrees{$Ped} }) {
	    my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
            my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
            if (($Mom != AttributeMissing) && ($Dad != AttributeMissing)) {
		$Seen{ $Mom . "+" . $Dad }++;
	    }
            if (($Dad == AttributeMissing) && ($Mom == AttributeMissing)) {
		$PedMsg{$Ped}{Fdrs}++;
		if ($Pedigrees{$Ped}{$Ind}{MkC} eq FALSE) { $PedMsg{$Ped}{FdrsMM}++; }
            } else {
		$PedMsg{$Ped}{NoFdrs}++;
	    }
            if (($Pedigrees{$Ped}{$Dad}{MkC} eq FALSE) && ($Pedigrees{$Ped}{$Mom}{MkC} eq FALSE)) {
                $PedMsg{$Ped}{MsgPPM}++;
            }
            if ($Pedigrees{$Ped}{$Ind}{MkC} eq FALSE) { $PedMsg{$Ped}{MsgMkr}++; }
            if ($Pedigrees{$Ped}{$Ind}{nPs} == AttributeMissing) { $PedMsg{$Ped}{MsgSib}++; }
            if ($Pedigrees{$Ped}{$Ind}{nMs} == AttributeMissing) { $PedMsg{$Ped}{MsgSib}++; }
            if ($Pedigrees{$Ped}{$Ind}{Aff} == AttributeMissing)    { $PedMsg{$Ped}{MsgAff}++; }
            if ($Pedigrees{$Ped}{$Ind}{Sex} == AttributeMissing) { $PedMsg{$Ped}{MsgSex}++; }
        }
        push @PedNuk, keys(%Seen) + 0;
        if ($PedMsg{$Ped}{Fdrs} == $PedMsg{$Ped}{FdrsMM}) {
	    print "Pedigree $Ped of " . $PedMsg{$Ped}{Fdrs} . " founders and " . $PedMsg{$Ped}{NoFdrs} . " non-founders has no founder markers!\n";
	}
    }

    @PedSiz = sort numerically @PedSiz;
    print "Number of pedigrees:" . @PedSiz . "\n";
    $Avg = sprintf("%.2f", (sum @PedSiz) / @PedSiz);
    @PedNuk = sort numerically @PedNuk;
    print "...Sizes: min:" . $PedSiz[0] . " max:" . $PedSiz[-1] . " med:" . $PedSiz[ $#PedSiz / 2 ] . " avg: $Avg\n";
    $Avg = sprintf("%.2f", (sum @PedNuk) / @PedNuk);
    print "...Nuclear family counts: min:" . $PedNuk[0] . " max:" . $PedNuk[-1] . " med:" . $PedNuk[ $#PedNuk / 2 ] . " avg: $Avg\n";

    for my $i (0 .. 7) {
	my @Counts = ();
	for my $Ped (keys %Pedigrees) {
	    push @Counts, $PedMsg{$Ped}[$i] + 0;
	}
	@Counts = sort numerically @Counts;
	print "..." . $CountsLabels{$i} . ": ";
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
	print "Family $Ped Depth: max:"
	    . $Ancestors[-1] . " med:"
	    . $Ancestors[ $#Ancestors / 2 ]
	    . " avg: $Avg\n";
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
# and both parents being treated as if unaffected, which will fall first into 
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

	# Adopt single case/control individuals
	my $memberCount = scalar(keys %{ $Pedigrees{$Ped} });
	if ($memberCount == 1) {
	    # &&& TBS
	}

        # Qualify the family for inclusion in trio buckets...
	if (($memberCount != 3) && ($memberCount != 4)) {
	    print "Copy $memberCount-individual pedigree $Ped directly\n";
	    next;
	}
	# Make a bucket prefix from the combined affection status of the parents
	# &&& TBS &&&

	my $skip = "";
	for my $Ind (sort keys %{ $Pedigrees{$Ped} }) {
	    if ((($Pedigrees{$Ped}{$Ind}{Dad} == AttributeMissing) && ($Pedigrees{$Ped}{$Ind}{Aff} != $Unaffected)) ||
		(($Pedigrees{$Ped}{$Ind}{Dad} != AttributeMissing) && ($Pedigrees{$Ped}{$Ind}{Aff} != $Affected))) {
#		$skip = "Ind $Ind has Dad ".$Pedigrees{$Ped}{$Ind}{Dad}." and affection status ".$Pedigrees{$Ped}{$Ind}{Aff};
		last;
	    }
	}
	if ($skip ne "") {
	    print "Copy bad affection pattern pedigree $Ped directly\n";
	    next;
	}

	# Get the family genotype bucket for each marker pair
	for my $i (0..$PairCount) {
	    my @bucketList = ();
	    # Get a trio bucket for each child in the family
	    for my $Ind (sort keys %{$Pedigrees{$Ped}}) {
		my $Dad = $Pedigrees{$Ped}{$Ind}{Dad};
		if ($Dad != AttributeMissing) {
		    my $ChildAlleles = $Pedigrees{$Ped}{$Ind}{Mks}[$i];
		    ($ChildAlleles eq ' 2 1') and $ChildAlleles = ' 1 2';
		    my $DadAlleles = $Pedigrees{$Ped}{$Dad}{Mks}[$i];
		    ($DadAlleles eq ' 2 1') and $DadAlleles = ' 1 2';
		    my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
		    my $MomAlleles = $Pedigrees{$Ped}{$Mom}{Mks}[$i];
		    ($MomAlleles eq ' 2 1') and $MomAlleles = ' 1 2';
#                        print "Get trio bucket for $MomAlleles $DadAlleles $ChildAlleles\n";
		    my $TrioBucket = $TrioBuckets{$MomAlleles}{$DadAlleles}{$ChildAlleles};
		    if ($TrioBucket eq "") {
			print "Couldn't find a bucket for pedigree $Ped, marker ".$Loci[$i+1].", M/D/C $MomAlleles/$DadAlleles/$ChildAlleles!\n";
			exit;
		    }
#			print "Pedigree $Ped / Marker ".$Loci[$i+1]." child $Ind (".$MomAlleles."-".$DadAlleles."-".$ChildAlleles.") gets bucket $TrioBucket\n";
		    push @bucketList, $TrioBucket;
		}
	    }
	    my $FullBucket = join ("-", sort (@bucketList) )."_".$Loci[$i+1];
#		print "Ped $Ped goes into full bucket $FullBucket for marker ".$Loci[$i+1]."\n";
	    if (!defined ($Buckets{$FullBucket})) {
#                    print "Defining bucket [$FullBucket]\n";
		$Buckets{$FullBucket} = 1;
	    } else {
#		    print "Bumping bucket [$FullBucket]\n";
		$Buckets{$FullBucket}++;
	    }
	}
    }
    foreach my $Bucket (sort keys (%Buckets)) {
	print "Bucket $Bucket has ".$Buckets{$Bucket}." entries\n";
    }
}

# Verify command line parameters
die "$#ARGV arguments supplied. $Usage" if ($#ARGV < 0);

my $ConfFile = shift;

loadConf($ConfFile);
my $companionFile = "datafile.dat";
if ($Directives{DF}[0] ne "") { $companionFile = $Directives{DF}[0]; }
loadCompanion($companionFile);
my $markersFile = "markers.dat";
if ($Directives{MK}[0] ne "") { $markersFile = $Directives{MK}[0]; }
loadMarkers($markersFile);
my $pedFile = "pedpost.dat";
if ($Directives{PD}[0] ne "") { $pedFile = $Directives{PD}[0]; }
assessPedigree($pedFile);

#loadPostPedigree($pedFile);
#bucketizePedigrees();

#checkRelations();
#perfStats();
simpleLoop();

exit;
