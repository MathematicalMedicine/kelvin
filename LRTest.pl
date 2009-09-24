#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#####################################
#
# $Id$
#
# by Bill Valentine-Cooper
#
# Copyright 2008, Nationwide Children's Research Institute.  All
# rights reserved.  Permission is hereby given to use this software
# for non-profit educational purposes only.
#
# Perform LR comparison tests on a kelvin dynamic grid configuration. Makes 
# extensive use of GNU baseutils, i.e. grep, sort, uniq. Expects TEST_KELVIN
# environment variable to point to version of kelvin to run. Requires kelvin
# V0.38 or later to support new directives and command-line directives.
#
# 1. Run kelvin with the kelvin.conf in the current directory, with the addition
# of the SurfaceFile directive.
#
# 2. Parse the SurfaceFile output header (starts with some spacing variation of 
# "# HLOD') to identify the trait space variables. All trait space variables must 
# have equivalent directives in %nameDirectives.
#
# 3. Find rows in the SurfaceFile output that have min and max values for each of
# the trait space variables. These will be used as test points.
#
# 4. Find 10? more randomly-chosen rows in the SurfaceFile output that are not the
# same as the ones already chose. These will be used as additional test points.
#
# 5. For each test point, run kelvin with the same kelvin.conf, with the addition 
# of the FixedGrids directive and all of the directives needed to choose that single 
# trait space test point.
#
# 6. Verify that all of the SurfaceFile outputs generated by the fixed grid test runs 
# are in the dynamic grid SurfaceFile output. If not, there's a problem.
#

my $EXTRACHOICES = 10; # Number of random choices excluding min/max choices
sub numerically { $a <=> $b } # Sort things numerically

# Associate column headers with configuration file directives
my %nameDirectives = (
		      'Theta' => 'Theta', 
		      'Alpha' => 'Alpha',
		      'DGF' => 'DiseaseGeneFrequency', 
		      'LC0PV-DD' => 'Penetrance DD',
		      'LC0PV-Dd' => 'Penetrance Dd',
		      'LC0PV-dD' => 'Penetrance dD',
		      'LC0PV-dd' => 'Penetrance dd',
		      'LC0DoFV-DD' => 'DegreesOfFreedom DD',
		      'LC0DoFV-Dd' => 'DegreesOfFreedom Dd',
		      'LC0DoFV-dD' => 'DegreesOfFreedom dD',
		      'LC0DoFV-dd' => 'DegreesOfFreedom dd',
		      'LC0MV-DD' => 'Mean DD',
		      'LC0MV-Dd' => 'Mean Dd',
		      'LC0MV-dD' => 'Mean dD',
		      'LC0MV-dd' => 'Mean dd',
		      'Dprime' => 'DPrime',
		      'SD' => 'StandardDev',
		      'Thresh' => 'Threshold',
		      'MkIdx' => 'Dummy', # Dummy directive (no action)
		      'PosIdx' => 'Dummy', # "
		      );

my  @headerNames = ();
sub parse_header {
    my @words = split;
    shift @words; # Ditch HLOD
    for my $word (@words) {
	if ($word =~ /(\S+)\((.*,.*)\)/) {
	    for my $partial (split(",",$2)) {
		die "Directive unknown for group column header ".$1."-".$partial."\n"
		    if (!defined($nameDirectives{$1."-".$partial}));
		push @headerNames, $1."-".$partial;
	    }
	} else {
	    die "Directive unknown for column header $word\n"
		if (!defined($nameDirectives{$word}));
	    push @headerNames, $word;
	}
    }
}

# Run kelvin on the config file with SurfaceFile directive
print "Generating dynamic-grid trait space baseline.\n";
my $commandLine='$TEST_KELVIN kelvin.conf --SurfaceFile LRTest.Dyn >&LRTest.Out';
(system($commandLine) == 0)
    or die "Couldn't run kelvin ($commandLine)\n";

# First pass thru file, pull-out min and max rows...
open IN,"LRTest.Dyn";

while (<IN>) {
    s/^\s*//g;       # Trim leading whitespace
    &parse_header($_) if (/\#\w*HLOD/);
    s/\s*\#.*//g;    # Trim comments
    next if (/^$/);  # Drop empty lines
    last;
}
my @words = split;
my $paramCnt = ($#words)-1;
my $HLOD = shift(@words); # Current HLOD
my @PsHLOD = (); # Each parameter's HLOD
my @PsTSV = (); # Each parameter's trait space values
my @PsLine = (); # Each parameter's original file line number
# Mins for each parameter are in element corresponding to parameter ordinality, and
# maxes are in element corresponding to parameter ordinality + number of parameters.
for my $i (0..$paramCnt) {
    $PsHLOD[$i] = $PsHLOD[$i+$paramCnt] = $HLOD;
    $PsTSV[$i] = $PsTSV[$i+$paramCnt] = [ @words ];
    $PsLine[$i] = $PsLine[$i+$paramCnt] = 0;
}
my $line = 0;
while (<IN>) {
    $line++;
    @words = split;
    $HLOD = shift(@words);
    # Pull the entry with the lowest/highest values for each trait space parameter
    for my $i (0..$paramCnt) {
	if ($words[$i] < $PsTSV[$i][$i]) {
	    $PsHLOD[$i] = $HLOD;
	    $PsTSV[$i] = [ @words ];
	    $PsLine[$i] = $line;
	}
	if ($words[$i] > $PsTSV[$i+$paramCnt][$i]) {
	    $PsHLOD[$i+$paramCnt] = $HLOD;
	    $PsTSV[$i+$paramCnt] = [ @words ];
	    $PsLine[$i+$paramCnt] = $line;
	}
    }
}
close IN;
my %avoidLines = ();
for my $i (0..$paramCnt) {
#    print "Parameter $i min at line ".$PsLine[$i]." of ".$PsTSV[$i][$i]." has HLOD ".
#	$PsHLOD[$i]."/TSV ".Dumper($PsTSV[$i])."\n";
    $avoidLines{$PsLine[$i]} = 1;
#    print "Parameter $i max at line ".$PsLine[$i+$paramCnt]." of ".
#	$PsTSV[$i+$paramCnt][$i]." has HLOD ".$PsHLOD[$i+$paramCnt].
#	"/TSV ".Dumper($PsTSV[$i+$paramCnt])."\n";
    $avoidLines{$PsLine[$i+$paramCnt]} = 1;
}
# Second pass thru file choosing N random ones not already chosen! This
# can be a fascinating problem if you're trying to find the 
# "cardinal intersection" of a set of columns. Unfortunately it would
# take time to figure out, and we're (always) in a hurry.
my @randLines = ();
my $choices = 0;
while ($choices < $EXTRACHOICES) {
    my $choice = int(rand($line));
    if (!defined($avoidLines{$choice})) {
	push @randLines, $choice;
	$choices++;
    }
}
@randLines = sort numerically @randLines;
#print "...and grab these lines too: ".Dumper(\@randLines)."\n";
open IN,"LRTest.Dyn";
$line = 0;
my $nextLine = shift(@randLines);
my $offset = $paramCnt * 2;
while (<IN>) {
    $line++;
    if ($line == $nextLine) {
	@words = split;
	$PsHLOD[$offset] = shift(@words);
	$PsTSV[$offset] = [ @words ];
	$PsLine[$offset] = $line;
	$offset++;
	$nextLine = shift(@randLines);
	last if (!defined($nextLine));
    }
}
close IN;

# Clean-out any earlier results. Allow this to fail.
system('rm LRTest-* >&/dev/null');

# Generate all of the kelvin configuration variations and run them
my $searchFormat = " " . ("%5.3f " x $paramCnt) . "%d";
for my $i (0..($offset - 1)) {
    my $commandLine = '$TEST_KELVIN kelvin.conf --FixedModels';
#    print "Test from line ".$PsLine[$i]." has HLOD ".
#	$PsHLOD[$i]."/TSV ".Dumper($PsTSV[$i])."\n";
    for my $j (0..$paramCnt) {
#	print "|".$nameDirectives{$headerNames[$j]}." ".$PsTSV[$i][$j]."\n";
	$commandLine .= " --".$nameDirectives{$headerNames[$j]}." ".$PsTSV[$i][$j];
    }
    $commandLine .= " --SurfaceFile LRTest-$i.Dat >& LRTest-$i.Out";
    $commandLine =~ s/--Dummy\s+[0-9]+/ /g if ($commandLine =~ /Dummy/); # Lose our dummy directives
    print "Generating fixed-grid trait space test point $i from dynamic grid output line ".$PsLine[$i].".\n";
#    print "Execute [$commandLine]\n";
    (system($commandLine) == 0) or die "Couldn't run \'$commandLine\'\n";
# Filter-out unrequested results, e.g. same parameters different marker, 0.5 Theta, etc.
    my $searchLine = sprintf($searchFormat, @{ $PsTSV[$i] }, $PsTSV[$i][$paramCnt]);
#    print "Limit by [$searchLine]\n";
    $commandLine = "grep \'$searchLine\' LRTest-$i.Dat >LRTest-$i.Fix";
    (system($commandLine) == 0) or die "Couldn't run \'$commandLine\'\n";
}

# Now concatenate, uniq-ify and compare all of the output, then let our driver do
# the differences and react accordingly.
(system('cat LRTest-*.Fix | sort | uniq > LRTest.Fix') == 0) or
    die "Cannot concatenate, sort or uniq all fixed-grid results\n";

# Find each fixed line in the dynamic output (using a possibly rounded HLOD).
print "Finding fixed results in dynamic results";
open IN,"LRTest.Fix";
while (<IN>) {
    chomp;
    print ".";
    if ($_ =~ /^([ \-][0-9]\.[0-9]{3})/) {
	my $old = $1;
# Heavy on the rounding slop. Go ahead and ask Sang for more surface precision!
	my $new = sprintf("[%5.3f |%5.3f |%5.3f |%5.3f |%5.3f |%5.3f |%5.3f |%5.3f |%5.3f |%5.3f ]",
			  $1-0.005, $1-0.004, $1-0.003, $1-0.002, $1-0.001, $1, $1+0.001, $1+0.002, $1+0.003, $1+0.004, $1+0.005);
# Deal with +/- zero
	$new =~ s/[ \-]0.000/ 0.000|-0.000/g;
#	print "HLOD is [$old], using $new\n";
	s/$old/$new/;
    }
    (system("grep \'".$_."\' LRTest.Dyn >&LRTest.grep") == 0) or
	die "Couldn't find line \'$_\' in LRTest.Dyn\n";
}
close IN;
print "done!\n";
