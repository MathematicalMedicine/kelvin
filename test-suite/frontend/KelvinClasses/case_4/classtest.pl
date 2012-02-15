#!perl -w
use strict;
use KelvinDataset;
use KelvinFamily;

my $dataset;
my $copy;
my $family;
my $copyfamily;
my $individual;

$dataset = KelvinDataset->new ({pedigreefile => 'ped.pre', locusfile => 'locus.dat',
				frequencyfile => 'freq.dat', mapfile => 'map.dat'})
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");
$dataset->misordered or exit (0);
$copy = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy->write ({pedigreefile => 'ped.new', locusfile => 'locus.new', mapfile => 'map.new',
	       freqfile => 'freq.new', backup => 0});

while ($family = $dataset->readFamily) {
    $copyfamily = $family->map ($copy);
    $copyfamily->write;
}
defined ($family) or die ("readFamily failed, $KelvinDataset::errstr\n");
