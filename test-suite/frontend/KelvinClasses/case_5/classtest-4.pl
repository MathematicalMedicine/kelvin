#!perl -w
use strict;
use KelvinDataset;

my $dataset;
my $copy;

$dataset = KelvinDataset->new ({locusfile => 'locus-4.dat', frequencyfile => 'freq-4.dat',
				mapfile => 'map-4.dat'})
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");
$dataset->misordered or die ("expected dataset to be misordered\n");
$copy = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy->write ({locusfile => 'locus.new', mapfile => 'map.new',
	       freqfile => 'freq.new', backup => 0});
exit (0);
