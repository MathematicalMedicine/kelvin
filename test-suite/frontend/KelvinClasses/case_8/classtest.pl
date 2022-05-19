#!perl -w
# Copyright (C) 2018, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
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
