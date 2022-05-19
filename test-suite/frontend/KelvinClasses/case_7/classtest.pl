#!perl -w
# Copyright (C) 2017, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
use strict;
use KelvinDataset;
use KelvinConfig;

# Test if we can use a KelvinConfig object to initialize a KelvinDataset, and
# also whether or not KelvinDataset::expandTraitPositions works as designed

my $config = KelvinConfig->new("test.conf")
        or die("KelvinConfig new failed: $KelvinConfig::errstr\n");
# assert new config worked
(${$config->isConfigured("PedigreeFile")}[0] eq "ped.post")
        or die("could not read config directive: $KelvinConfig::errstr\n");

my $dataset = KelvinDataset->new($config)
        or die("KelvinDataset new failed: $KelvinDataset::errstr\n");
# assert dataset exists and can read from our datafiles
(defined($dataset->readIndividual()))
        or die("could not readIndividual: $KelvinDataset::errstr\n");

my $traitpositions = $dataset->expand_traitpositions(
        ${$config->isConfigured("TraitPositions")}[0]);

# assert each position is present
my $basepositions = [0, 2, 4, 5, 100, 101, 102, 103, 104, 105, 106, 107, 108,
        109];
for (my $i = 0; $i < scalar(@$basepositions); $i++) {
    ($$traitpositions[$i] == $$basepositions[$i])
            or die("expanded TraitPositions mismatch at index $i");
}
