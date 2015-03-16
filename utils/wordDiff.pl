#!/usr/bin/env perl

# FIXME: why isn't 'use strict;' present?
use warnings;

#
# Word-oriented file comparison with case folding and a numeric slop factor.
#
# Bill Valentine-Cooper - William.Valentine-Cooper@NationwideChildrens.org
#
# Copyright 2010, The Research Institute at Nationwide Children's Hospital
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
# verify that files exist before we proceed
unless (-e $file_LEFT) {
    print "File $file_LEFT not found.\n";
    exit (1);
}
unless (-e $file_RIGHT) {
    print "File $file_RIGHT not found.\n";
    exit (1);
}
open LEFT, $file_LEFT; open RIGHT, $file_RIGHT;
$whole_LEFT = do { local $/; <LEFT> }; # Briefly make the special variable for line delimiter undefined...
$whole_RIGHT = do { local $/; <RIGHT> }; # ...so we can slurp-in the whole file at once.
close LEFT; close RIGHT;
$whole_LEFT =~ s/Version V.*//g; # Ignore comparisons of Version strings
$whole_RIGHT =~ s/Version V.*//g; # Ignore comparisons of Version strings
$whole_LEFT =~ s/\$Id.* \$//g; # Ignore comparisons of SVN or CVS version Ids
$whole_RIGHT =~ s/\$Id.* \$//g; # Ignore comparisons of SVN or CVS version Ids
# Split the files into lists delimited by whitespace, commas, newlines and parens, and keep delimiters
@chunks_LEFT = split (/([\s\n]+|[,\(\)])/, $whole_LEFT);
@chunks_RIGHT = split (/([\s\n]+|[,\(\)])/, $whole_RIGHT);
# Compare all of the tokens and delimiters
for ($i = 0; $i < $#chunks_LEFT; $i++) {
    $left = $chunks_LEFT[$i]; $right = $chunks_RIGHT[$i];
    next if ($left =~ /V\S+/ && $right =~ /V\S+/); # Ignore comparisons of version numbers (e.g. V2.3)
    next if ($left =~ /^\s*$/ && $right =~ /^\s*$/); # chunks of contiguous whitespace always match
    if ($left ne $right) {
	# Not a textual match
	if ($left =~ /^[-+]?[0-9]+[\.]?[0-9]*([Ee][+-][0-9]*)?$/) {
	    # It's numeric, try simple equality to avoid e+/-00 issues.
	    if ($left != $right) {
		# Not simply equal, get the deviation and compare to maximum allowable
		if ($left != 0) {
		    $offByFactor = ($left-$right)/$left;
		} else {
		    print "Zero value prevents proper comparison\n";
		    exit (1);
		}
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
