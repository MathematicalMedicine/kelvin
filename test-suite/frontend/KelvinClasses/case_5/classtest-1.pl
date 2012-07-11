#!perl -w
use strict;
use KelvinDataset;

my $dataset;
my $copy;

$dataset = KelvinDataset->new ({locusfile => 'locus-1.dat', frequencyfile => 'freq-1.dat'})
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");
$dataset->misordered or die ("expected dataset to be misordered\n");
$copy = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy->write ({locusfile => 'locus.new', freqfile => 'freq.new', backup => 0});
exit (0);
