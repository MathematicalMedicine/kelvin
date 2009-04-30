#!/usr/bin/perl -w
#
# Word-oriented file comparison with case folding and a numeric slop factor.
#
# Bill Valentine-Cooper - William.Valentine-Cooper@NationwideChildrens.org
#
# Copyright 2008, The Research Institute at Nationwide Children's Hospital
# All rights reserved. Permission is granted to use this software for
# non-profit educational purposes only.
#

if ($#ARGV < 2) {
    print "Usage: ".$ARGV[0]." <file> <file> <maximum numeric error factor, 0 for equality>\n";
    exit (2);
}
# Read in both files into scalars in their entirety so we can use split on them.
$file_LEFT = shift;
$file_RIGHT = shift;
$maxOffByFactor = shift;
open LEFT, $file_LEFT; open RIGHT, $file_RIGHT;
$whole_LEFT = do { local $/; <LEFT> }; # Briefly change newline to NULL
$whole_RIGHT = do { local $/; <RIGHT> }; # Briefly change newline to NULL
close LEFT; close RIGHT;
# Split the files into lists delimited by whitespace commas, newlines and parens, and keep delimiters
@chunks_LEFT = split (/([\s,\n\(\)])/, $whole_LEFT);
@chunks_RIGHT = split (/([\s,\n\(\)])/, $whole_RIGHT);
# Compare all of the tokens and delimiters
for ($i = 0; $i < $#chunks_LEFT; $i++) {
    $left = $chunks_LEFT[$i]; $right = $chunks_RIGHT[$i];
    if ($left ne $right) {
	# Not a textual match
	if ($left =~ /^[-+]?[0-9]+[\.]?[0-9]*([Ee][+-][0-9]*)?$/) {
	    # It's numeric, try simple equality to avoid e+/-00 issues.
	    if ($left != $right) {
		# Not simply equal, get the deviation and compare to maximum allowable
		$offByFactor = ($left-$right)/$left;
#		print "Off by $offByFactor vs limit of $maxOffByFactor\n";
		if (abs($offByFactor) > $maxOffByFactor) {
		    print "Difference of factor of $offByFactor between $left and $right ".
			"is over limit of $maxOffByFactor\n";
		    exit (1);
		}
	    }
	} else {
	    # Non-numeric match problem, try case folding
	    if (lc($left) ne lc($right)) {
		print "Difference between [$left] and [$right]\n";
		exit (1);
	    }
	}
    }
}
exit 0;
