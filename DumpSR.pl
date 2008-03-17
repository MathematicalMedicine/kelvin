#!/usr/bin/perl
use lib "./";
use Tpl;

$fileName = shift();

my $traitTPLFormat = "sf#f#f#f#f#f#";
my $traitFileFormat = "trait.tpl";
my $markerTPLFormat = "siiA(s)f";
my $markerFileFormat = "marker.tpl";
my $alternativeTPLFormat = "siff#f#f#f#f#f#";
my $alternativeFileFormat = "alternative.tpl";
my $pedigree, $chromosome, $traitPosition, @lDT0s, @lDT1s, @lDT2s, @lDT3s, @lDT4s, @lDT5s;

if ($fileName =~ /$alternativeFileFormat/) {
    dumpAlternative();
} else {
    if ($fileName =~ /$markerFileFormat/) {
	dumpMarker();
    } else {
	if ($fileName =~ /$traitFileFormat/) {
	    dumpTrait();
	}
    }
}

sub dumpMarker() {
    $tpl = Tpl->tpl_map($markerTPLFormat, \$pedigree, \$chromosome, \$markerCount, \$markerName,
			\$likelihood);
    $tpl->tpl_load($fileName);
    $tpl->tpl_unpack(0);
    $value = sprintf("%.6g", $likelihood);
    print "Pedigree $pedigree, chromosome $chromosome, $markerCount markers, likelihood $value\n";
    for ($i=0;$i<$markerCount;$i++) {
	$tpl->tpl_unpack(1);
	print "$i: $markerName\n";
    }
}

sub dumpTrait() {
    my @lDT;
    $lDT[0][0] = 0;
    $tpl = Tpl->tpl_map($traitTPLFormat, \$pedigree,
			\@{$lDTs[0]}, 275, \@{$lDTs[1]}, 275, \@{$lDTs[2]}, 275, \@{$lDTs[3]}, 275,
			\@{$lDTs[4]}, 275, \@{$lDTs[5]}, 275);
    $tpl->tpl_load($fileName);
    $tpl->tpl_unpack(0);
    print "Pedigree $pedigree\n";
    for ($i=0;$i<6;$i++) {
	for ($j=0;$j<275;$j++) {
	    if (($j % 4) == 0) {
		print "\n";
	    }
	    $value = sprintf("%.6g", $lDTs[$i][$j]);
	    print "($i,$j): $value\t";
	}
    }
    print "\n";
}

sub dumpAlternative() {
    $tpl = Tpl->tpl_map($alternativeTPLFormat, \$pedigree, \$chromosome, \$traitPosition,
			\@{$lDTs[0]}, 275, \@{$lDTs[1]}, 275, \@{$lDTs[2]}, 275, \@{$lDTs[3]}, 275,
			\@{$lDTs[4]}, 275, \@{$lDTs[5]}, 275);
    $tpl->tpl_load($fileName);
    $tpl->tpl_unpack(0);

    print "Pedigree $pedigree, chromosome $chromosome, position $traitPosition\n";
    for ($i=0;$i<6;$i++) {
	for ($j=0;$j<275;$j++) {
	    if (($j % 4) == 0) {
		print "\n";
	    }
	    $value = sprintf("%.6g", $lDTs[$i][$j]);
	    print "($i,$j): $value\t";
	}
    }
    print "\n";
}
