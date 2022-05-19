# Copyright (C) 2008, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
open IN, shift();
$name = shift();
open DATA, ">$name.dat";
open CMD, ">$name.cmd";
print CMD "set terminal gif size 1024,768;\n";
print CMD "set output \"$name.gif\";\n"; 
print CMD "set title \"Computational Complexity for $name\"\n";
print CMD "set xlabel \"Chromosome Position (cM)\"\n";
print CMD "set ylabel \"Pedigree (in sequence)\"\n";
print CMD "set pm3d map;\n";
%Uniques = {}; %Similars = {};
$lowPos = 1000; $highPos = 0;
while (<IN>) {
    if (/Ped ([0-9]+)\([0-9]+\) has/) {
	$pedigree = $1;
#	print "Pedigree $1 has ";
    } else {
	if (/\b*Ped has total ([0-9]*) unique pp groups\, ([0-9]*)/) {
	    $Uniques{$pedigree} = $1; $Similars{$pedigree} = $2;
#	    print "$1 unique and $2 similar\n";
	} else {
	    if (/POS ([0-9]*\.?[0-9]*) has .*total ([0-9]*)./) {
		$position = $1; $total = $2;
		if ($position < $lowPos) {
		    $lowPos = $position;
		}
		if ($position > $highPos) {
		    $highPos = $position;
		}
#		print "Grand total for position $1 is $2\n";
		$i = 0;
		for $pedigree (sort(keys %Uniques)) {
		    if ($pedigree > 0) {
			print DATA $position." ".(++$i)." ".($Uniques{$pedigree}+$Similars{$pedigree})."\n";
			if (($Uniques{$pedigree}+$Similars{$pedigree}) > 1e+9) {
			    print CMD "set label $i \"$pedigree\" at -20,$i,1\n";
			}
		    }
		}
		print DATA "\n";
	    }
	}
    }
}
close IN;
close DATA;
print CMD "set xrange [$lowPos:$highPos];\n";
print CMD "set yrange [1:".((scalar keys %Uniques)-1)."];\n";
#print CMD "set zrange [0:2e+9];\n";
print CMD "splot \"$name.dat\";\n";
close CMD;
system("gnuplot $name.cmd");
