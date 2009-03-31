# Digits then whitespace then string then string then signed decimal is pretty convincing, so
# take the rest of the signed decimals on faith so we don't go crazy with the expression.
# Chr Trait=1 Marker=2 Position=3 MOD=4 D11=5 Theta(M=6,F) Alpha=7 DGF=8 LC0PV(DD=9,Dd=10,dd=11)
if (/^\d+\s+(\S+)\s+(\S+)\s+([-+]?[0-9]*\.?[0-9]+)\s+(\S+)\s+(\S+)\s+\((\S+),\S+\)\s+(\S+)\s+(\S+)\s+\((\S+),(\S+),(\S+)\)/) {
    open OUT,">single-$2.conf";
    print OUT "MD mod-$2.out\nLD $5\nTh $6\nAL $7\nGF $8\nDD $9\nDd $10\ndd $11\n";
    close OUT;
    system($ENV{TEST_KELVIN} . " single-$2.conf >single-$2.out 2>&1") and die "Failed to run kelvin on marker $2\n";
    system("grep $2 mod.out >dynamic-$2.out; grep $2 mod-$2.out >fixed-$2.out") and die "Failed to extract results for marker $2\n";
    system("diff dynamic-$2.out fixed-$2.out") and die "Dynamic and fixed differ for marker $2\n";
}
