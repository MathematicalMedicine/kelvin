#!/usr/bin/env perl

$StudyId = shift();
$ChromosomeNo = shift();
$MarkerCount = shift();

open SQL,">SMRT_Inserts.sql";

foreach $file (glob("kelvin.out-*")) {
    open IN,$file;
    $file =~ /.*-(\d+)/;
    $ped = $1;
    $prevPos = "";
    $prevSec = 0;
    while (<IN>) {
	if (/\s+\@((\d+)d)?((\d+)h)?((\d+)m)?(\d+)s.*trait locus at (\d+\.?\d*) /) {
	    if ($5 eq "") {
		$pos = $8; $sec = $7;
	    } elsif ($4 eq "") {
		$pos = $8;  $sec = $7+($6*60.0);
	    } elsif ($2 eq "") {
		$pos = $8; $sec = $7+($6*60.0)+($4*3600.0);
	    } else {
		$pos = $8; $sec = $7+($6*60.0)+($4*3600.0)+($2*3600*24);
	    }
	    if ($prevPos ne "") {
		print SQL "Insert into SingleModelRuntimes (StudyId, PedigreeSId, ChromosomeNo, PedTraitPosCM, MarkerCount, SingleModelRuntime) values ".
		    "($StudyId, \'$ped\', $ChromosomeNo, $prevPos, $MarkerCount, ".($sec - $prevSec).");\n";
	    }
	    $prevPos = $pos;
	    $prevSec = $sec;
	}
	if (/\s+\@((\d+)d)?((\d+)h)?((\d+)m)?(\d+)s.*Analysis complete/) {
	    if ($5 eq "") {
		$sec = $7;
	    } elsif ($4 eq "") {
		$sec = $7+($6*60.0);
	    } elsif ($2 eq "") {
		$sec = $7+($6*60.0)+($4*3600.0);
	    } else {
		$sec = $7+($6*60.0)+($4*3600.0)+($2*3600*24);
	    }
	    print SQL "Insert into SingleModelRuntimes (StudyId, PedigreeSId, ChromosomeNo, PedTraitPosCM, MarkerCount, SingleModelRuntime) values ".
		"($StudyId, \'$ped\', $ChromosomeNo, $prevPos, $MarkerCount, ".($sec - $prevSec).");\n";
	}
    }
    close IN;
    $pedOrd++;
}

close SQL;
