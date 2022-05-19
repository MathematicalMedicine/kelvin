# Copyright (C) 2009, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
# 
# Digits then whitespace then string then string then signed decimal is pretty convincing, so
# take the rest of the signed decimals on faith so we don't go crazy with the expression.
# Chr Trait=1 Marker=2 Position=3 MOD=4 D11=5 Theta(M=6,F) Alpha=7 DGF=8 LC0PV(DD=9,Dd=10,dd=11)
if (/^\d+\s+(\S+)\s+(\S+)\s+([-+]?[0-9]*\.?[0-9]+)\s+(\S+)\s+(\S+)\s+\((\S+),\S+\)\s+(\S+)\s+(\S+)\s+\((\S+),(\S+),(\S+)\)/) {
    open OUT,">single-$2.conf";
    print OUT "MODFile mod-$2.out\nLD $5\nTheta $6\nAlpha $7\nDiseaseGeneFrequency $8\nPenetrance DD $9\nPenetrance Dd $10\nPenetrance dd $11\n";
    close OUT;
    system($ENV{TEST_KELVIN} . " single-$2.conf >single-$2.out 2>&1") and die "Failed to run kelvin on marker $2\n";
    system("grep $2 mod.out >dynamic-$2.out; grep $2 mod-$2.out >fixed-$2.out") and die "Failed to extract results for marker $2\n";
    system("diff dynamic-$2.out fixed-$2.out") and die "Dynamic and fixed differ for marker $2\n";
}
