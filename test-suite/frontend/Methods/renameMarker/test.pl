#!/usr/bin/env perl
use warnings;
use KelvinDataset;

my $set;

# Case 1: frequency file only
$set = KelvinDataset->new ({frequencyfile => 'test1.freq'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker1", "newmarker1")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({frequencyfile => "case1.freq"})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");
$set = undef;


# Case 2: matching frequency and map files
$set = KelvinDataset->new ({frequencyfile => 'test1.freq', mapfile => 'test1.map'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker1", "newmarker1")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({frequencyfile => 'case2.freq', mapfile => 'case2.map'})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");
$set = undef;


# Case 3: map file with subset locus file 
$set = KelvinDataset->new ({locusfile => 'test2.locus', mapfile => 'test1.map'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker1", "newmarker1")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({locusfile => 'case3.locus', mapfile => 'case3.map'})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");
$set = undef;


# Case 4: matching map and frequency files with subset locus file
$set = KelvinDataset->new ({locusfile => 'test2.locus', mapfile => 'test1.map',
			    frequencyfile => 'test1.freq'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker1", "newmarker1")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({locusfile => 'case4.locus', mapfile => 'case4.map', frequencyfile => 'case4.freq'})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");
$set = undef;


# Case 5: locus file only
$set = KelvinDataset->new ({locusfile => 'test2.locus'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker2", "newmarker2")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({locusfile => 'case5.locus'})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");
$set = undef;


# Case 6: matching map and locus files
$set = KelvinDataset->new ({locusfile => 'test2.locus', mapfile => 'test2.map'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker2", "newmarker2")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({locusfile => 'case6.locus', mapfile => 'case6.map'})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");
$set = undef;


# Case 6: matching map, frequency and locus files
$set = KelvinDataset->new ({locusfile => 'test2.locus', mapfile => 'test2.map',
			frequencyfile => 'test2.freq'})
    or die ("KelvinDataset constructor failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker2", "newmarker2")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker8", "newmarker8")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->renameMarker ("marker20", "newmarker20")
    or die ("KelvinDataset renameMarker failed, $KelvinDataset::errstr\n");
$set->write ({locusfile => 'case7.locus', mapfile => 'case7.map', frequencyfile => 'case7.freq'})
    or die ("KelvinDataset write failed, $KelvinDataset::errstr\n");

exit (0);
