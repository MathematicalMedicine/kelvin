#!perl -w
use strict;
use KelvinDataset v1.9.0;

my $dataset;
my $ind;

$dataset = KelvinDataset->new ({pedigreefile => 'test.pre', locusfile => 'test.locus' })
    or die ("KelvinDataset new failed: $KelvinDataset::errstr\n");

open (REP, ">", "report.txt");
print (REP "#FID IID Genos HetGenos\n");
while ($ind = $dataset->readIndividual) {
    print (REP join (" ", $ind->pedid, $ind->indid, $ind->genotyped, $ind->heterozygous), "\n");
}
defined ($ind) or die ("readIndividual failed, $KelvinDataset::errstr\n");
close (REP);
exit (0);
