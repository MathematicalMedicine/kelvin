#!perl -w
# Copyright (C) 2012, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
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
