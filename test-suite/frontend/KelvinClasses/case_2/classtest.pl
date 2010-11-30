#!perl -w
use strict;
use KelvinDataset;
use KelvinFamily;

my $dataset;
my $copy1;
my $copy2;
my $family;
my $copyfamily;
my $individual;


$dataset = KelvinDataset->new ({pedigreefile => 'ped.pre', locusfile => 'locus.dat',
				frequencyfile => 'freq.dat', mapfile => 'map.dat'})
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");

$copy1 = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy1->writePedigreefile ({pedigreefile => 'ped.post.new.1', backup => 0})
    or die ("KelvinDataset writePedgreefile failed: $KelvinDataset::errstr\n");

$copy2 = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy2->writePedigreefile ({pedigreefile => 'ped.pre.new.1', backup => 0, premakeped => 1})
    or die ("KelvinDataset writePedgreefile failed: $KelvinDataset::errstr\n");

while ($family = $dataset->readFamily) {
    $copyfamily = $family->map ($copy1);
    $copyfamily->write;
    $copyfamily = $family->map ($copy2);
    $copyfamily->write;
}
defined ($family) or die ("readFamily failed, $KelvinDataset::errstr\n");


$dataset = KelvinDataset->new ({pedigreefile => 'ped.post', locusfile => 'locus.dat',
				frequencyfile => 'freq.dat', mapfile => 'map.dat'})
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");

$copy1 = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy1->writePedigreefile ({pedigreefile => 'ped.post.new.2', backup => 0})
    or die ("KelvinDataset writePedgreefile failed: $KelvinDataset::errstr\n");

$copy2 = $dataset->copy
    or die ("KelvinDataset copy failed, $KelvinDataset::errstr\n");
$copy2->writePedigreefile ({pedigreefile => 'ped.pre.new.2', backup => 0, premakeped => 1})
    or die ("KelvinDataset writePedgreefile failed: $KelvinDataset::errstr\n");

while ($family = $dataset->readFamily) {
    $copyfamily = $family->map ($copy1);
    $copyfamily->write;
    $copyfamily = $family->map ($copy2);
    $copyfamily->write;
}
defined ($family) or die ("readFamily failed, $KelvinDataset::errstr\n");
