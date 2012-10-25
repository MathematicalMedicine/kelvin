#!perl -w
use strict;
use KelvinDataset;
use KelvinFamily;

my $dataset;
my $copy;
my $family;
my $copyfamily;
my $individual;

$dataset = KelvinDataset->new ({pedigreefile => 'ped.post', locusfile => 'locus.dat',
				frequencyfile => 'freq.dat', mapfile => 'map.dat'})
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");
$copy = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy->deleteTrait ("Affection")
    or die ("KelvinDataset deleteTrait failed, $KelvinDataset::errstr\n");
$copy->addMarker ("Loci_5", {chr => 40, name => 'Loci_5', avgpos => 0.01},
                  {51 => 0.445057, 52 => 0.554943})
    or die ("KelvinDataset addMarker failed, $KelvinDataset::errstr\n");
$copy->write ({pedigreefile => 'ped.new', locusfile => 'locus.new', mapfile => 'map.new',
	       freqfile => 'freq.new', backup => 0});

while ($family = $dataset->readFamily) {
    $copyfamily = $family->map ($copy);
    $copyfamily->write;
}
defined ($family) or die ("readFamily failed, $KelvinDataset::errstr\n");
