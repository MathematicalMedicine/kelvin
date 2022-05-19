#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions;
use Data::Dumper;

$|=1; # Immediate output

# kdiff --help for documentation.

# Used globally
my $svn_version='$Id$';


my $are_different = 0;
my $use_found_files = 0;

# Argument- and option-related
my $verbose = 0;
my $help = 0;
my $pedfiles = 0; my $dup_pedfiles = 0;
my $freqfiles = 0;
my $mapfiles = 0;
my $path1 = ''; my $path2 = '';
my $pedfile1 = ''; my $pedfile2 = '';
my $locusfile1 = ''; my $locusfile2 = '';
my $freqfile1 = ''; my $freqfile2 = '';
my $mapfile1 = ''; my $mapfile2 = '';
my $configfile1 = './kelvin.conf'; my $configfile2 = './kelvin.conf';
my $freqFuzz = 0; # Might want to allow up to 0.02;
my $posFuzz = 0; # Might want to allow up to 0.001;
my $tabular = 0;

# Hash references to keep KelvinDataset happy
my ($data1ref, $data2ref);

# Primary Kelvin* workhorses
my ($config1, $config2);
my ($dataset1, $dataset2);

# Something I'd expected to find in Kelvin* but wasn't there...
my (%families1, %families2);

my $KELVIN_ROOT='no_kelvin_root';

GetOptions (
    'verbose' => \$verbose,
    'help|?' => \$help,
    'pedfiles' => \$pedfiles,
    'freqfiles' => \$freqfiles,
    'mapfiles' => \$mapfiles,
    'configfile1|c1=s' => \$configfile1,
    'configfile2|c2=s' => \$configfile2,
    'path1|p1=s' => \$path1,
    'path2|p2=s' => \$path2,
    'freq-fuzz=s' => \$freqFuzz,
    'pos-fuzz=s' => \$posFuzz,
    'tabular' => \$tabular,
    ) or pod2usage(2);

pod2usage(-noperldoc => 1, -verbose => 2) if ($help); # All source is shown if -noperdoc isn't specified!

pod2usage(-verbose => 0) if ($#ARGV > -1);

check_paths ();

print "Processing configuration files 1:\"$configfile1\" and 2:\"$configfile2\"\n" if $verbose;

check_paths ();

$use_found_files = 1 if (($pedfiles eq 0) and ($freqfiles eq 0) and ($mapfiles eq 0));
print "Processing all \"found\" data files\n" if ($verbose and $use_found_files);

# First we load (and thereby validate) everything

($config1 = KelvinConfig->new ($configfile1))
    or error ("Processing of configuration file 1:\"$configfile1\" failed: $KelvinConfig::errstr");
#print "Config 1 is ".Dumper($config1)."\n";
($config2 = KelvinConfig->new ($configfile2))
    or error ("Processing of configuration file 2:\"$configfile2\" failed: $KelvinConfig::errstr");
#print "Config 2 is ".Dumper($config2)."\n";

if ($use_found_files or $freqfiles) {
    # Construct path and name for frequency files using default since KelvinConfig doesn't.
    my $filename = defined ($config1->isConfigured ("FrequencyFile")) ?
	$ {$config1->isConfigured ("FrequencyFile")}[0] : "markers.dat";

    $filename = catfile($path1, $filename)
	if ((!file_name_is_absolute($filename)) && $path1);
    # Admitted strangeness here...I'm letting KelvinConfig give them the bad news.
    $$data1ref{FrequencyFile} = $filename if ((-e $filename) or $freqfiles);
	
    $filename = defined ($config2->isConfigured ("FrequencyFile")) ?
	$ {$config2->isConfigured ("FrequencyFile")}[0] : "markers.dat";

    $filename = catfile($path2, $filename)
	if ((!file_name_is_absolute($filename)) && $path2);
    $$data2ref{FrequencyFile} = $filename if ((-e $filename) or $freqfiles);

    $freqfiles = 1 if ($$data1ref{FrequencyFile} and $$data2ref{FrequencyFile});
}

if ($use_found_files or $mapfiles) {
    # Construct path and name for map files
    my $filename = $ {$config1->isConfigured ("MapFile")}[0];
    $filename = catfile($path1, $filename)
	if ((!file_name_is_absolute($filename)) && $path1);
    $$data1ref{MapFile} = $filename if ((-e $filename) or $mapfiles);
	
    $filename = $ {$config2->isConfigured ("MapFile")}[0];
    $filename = catfile($path2, $filename)
	if ((!file_name_is_absolute($filename)) && $path2);
    $$data2ref{MapFile} = $filename if ((-e $filename) or $mapfiles);

    $mapfiles = 1 if ($$data1ref{MapFile}) and ($$data2ref{MapFile});
}

if ($use_found_files or $pedfiles) {
    # Construct path and name for pedigree and locus files
    my $pedname = $ {$config1->isConfigured ("PedigreeFile")}[0];
    my $locusname = $ {$config1->isConfigured ("LocusFile")}[0];
    $pedname = catfile($path1, $pedname)
	if ((!file_name_is_absolute($pedname)) && $path1);
    $locusname = catfile($path1, $locusname)
	if ((!file_name_is_absolute($locusname)) && $path1);
    if (((-e $pedname) and (-e $locusname)) or $pedfiles) {
	$$data1ref{PedigreeFile} = $pedname;
	$$data1ref{LocusFile} = $locusname;
    }
    $pedname = $ {$config2->isConfigured ("PedigreeFile")}[0];
    $locusname = $ {$config2->isConfigured ("LocusFile")}[0];
    $pedname = catfile($path2, $pedname)
	if ((!file_name_is_absolute($pedname)) && $path2);
    $locusname = catfile($path2, $locusname)
	if ((!file_name_is_absolute($locusname)) && $path2);
    if (((-e $pedname) and (-e $locusname)) or $pedfiles) {
	$$data2ref{PedigreeFile} = $pedname;
	$$data2ref{LocusFile} = $locusname;
    }

    # Turn on option if we find them, but don't turn it off if we don't so we can error-out as needed.
    if (($$data1ref{PedigreeFile}) and ($$data1ref{LocusFile}) and
	($$data2ref{PedigreeFile}) and ($$data2ref{LocusFile})) {
	$pedfiles = 1;
	# If all files are the same, don't bother validating again. We'll still compare just to verify this code.
	$dup_pedfiles = 1 if ((($$data1ref{PedigreeFile}) eq ($$data2ref{PedigreeFile})) and
			      (($$data1ref{LocusFile}) eq ($$data2ref{LocusFile})));
    }
}

# Read and validate everything we have (except the pedfile contents)
if ($verbose or $use_found_files) {
    print "1: Validating".($use_found_files ? " found" : "")." files: ";
    for my $key (sort %{$data1ref}) {print "$key->".$$data1ref{$key}." " if (defined($$data1ref{$key})); }
    print "\n";
}
$dataset1 = KelvinDataset->new ($data1ref)
    or error ("1: Validation of referenced or defaulted data files failed: $KelvinDataset::errstr");
if ($verbose or $use_found_files) {
    print "2: Validating".($use_found_files ? " found" : "")." files: ";
    for my $key (sort %{$data2ref}) {print "$key->".$$data2ref{$key}." " if (defined($$data2ref{$key})); }
    print "\n";
}
$dataset2 = KelvinDataset->new ($data2ref)
    or error ("2: Validation of referenced or defaulted data files failed: $KelvinDataset::errstr");
#print "Dataset2 is ".Dumper($dataset2)."\n";

if ($pedfiles) {
    # Read, validate and describe the pedigrees
    my $family;
    my $totalInds = 0;
    while ($family = $dataset1->readFamily) {
	$$dataset1{origfmt} = $$family{origfmt};
	$totalInds += $$family{count};
	print "1: ".$family->pedtype." family ".$$family{pedid}." of ".$$family{count}." (".$$family{founders}."f/".$$family{nonfounders}."nf)\n" if $verbose;
	$families1{$$family{pedid}} = $family;
    }
    (defined ($family))
	or error ("1: Read of pedigree file \"".$$data1ref{PedigreeFile}."\" failed, $KelvinDataset::errstr");
    print "1: ".$$dataset1{origfmt}."-makeped format file with ".scalar(@{$$dataset1{markerorder}})." markers, ".
	scalar(keys %families1)." families and $totalInds individuals.\n" if $verbose;

    if ($dup_pedfiles) {
	$dataset2 = $dataset1;
	%families2 = %families1; # Don't really want to double storage, so fix this.
    } else {
	$totalInds = 0;
	while ($family = $dataset2->readFamily) {
	    # print "Family structure is ".Dumper($family)."\n";
	    $$dataset2{origfmt} = $$family{origfmt};
	    $totalInds += $$family{count};
	    print "2: ".$family->pedtype." family ".$$family{pedid}." of ".$$family{count}." (".$$family{founders}."f/".$$family{nonfounders}."nf)\n" if $verbose;
	    $families2{$$family{pedid}} = $family;
	}
	(defined ($family))
	    or error ("2: Read of pedigree in file \"".$$data2ref{PedigreeFile}."\" failed, $KelvinDataset::errstr");
	print "2: ".$$dataset2{origfmt}."-makeped format file with ".scalar(@{$$dataset2{markerorder}})." markers, ".
	    scalar(keys %families2)." families and $totalInds individuals.\n" if $verbose;
    }
}

if ($mapfiles) {
    # Describe and compare the map files
    print "1: ".$$dataset1{mapfunction}." map of ".scalar(@{$$dataset1{maporder}})." markers for chromosome ".$$dataset1{chromosome}." providing ".join(", ",@{$$dataset1{mapfields}})."\n" if $verbose;
    print "2: ".$$dataset2{mapfunction}." map of ".scalar(@{$$dataset2{maporder}})." markers for chromosome ".$$dataset2{chromosome}." providing ".join(", ",@{$$dataset2{mapfields}})."\n" if $verbose;

    if ($tabular) {
	print "#A,Path1,Path2,Marker,Attribute,Value1,Value2,Difference\n";
	print "#R,Path1,Path2,Marker,PreviousMarker,Attribute,Value1,PreviousValue1,Value2,PreviousValue2,Difference1,Difference2\n";
    }

    my %fieldhash1 = map { $_ => 1 } @{$$dataset1{mapfields}};
    my %fieldhash2 = map { $_ => 1 } @{$$dataset2{mapfields}};
    my @common_fields = ();
    for my $field (uniqua ((@{$$dataset1{mapfields}}, @{$$dataset2{mapfields}}))) {
	if (!defined($fieldhash1{$field})) {
	    print "1: Map doesn't provide $field, will not be compared!\n";
	    $are_different += 1;
	    next;
	}
	if (!defined($fieldhash2{$field})) {
	    print "2: Map doesn't provide $field, will not be compared!\n";
	    $are_different += 1;
	    next;
	}
	push @common_fields, $field;
    }

# The following won't work until 5.10 shows up on the cluster! And upon reflection, order doesn't matter (to us).
#    if (!(@{$$dataset1{maporder}} ~~ @{$$dataset2{maporder}})) {
#	print "2: Order of markers in maps is different!\n";
#	$are_different += 1;
#    }

    my %map1 = %{$$dataset1{markers}}; my %map2 = %{$$dataset2{markers}};
    # Compare markers by looping over the superset of keys (marker names) in maporder. This is safe because the map must define a superset of markers
    my $lastName = "";
    for my $name (uniqua ((@{$$dataset1{maporder}}, @{$$dataset2{maporder}}))) {
	if (!defined($map1{$name})) {
	    print "1: Marker \"$name\" not found in map, skipping!\n";
	    $are_different += 1;
	    next;
	}
	if (!defined($map2{$name})) {
	    print "2: Marker \"$name\" not found in map, skipping!\n";
	    $are_different += 1;
	    next;
	}
    }
    # Now we have to go over the common markers in map order, so either list will work provided we ditch entries not present in the other
    for my $name (@{$$dataset1{maporder}}) {
	next if (!defined($map2{$name}));
	for my $field (@common_fields) {
	    if (($field =~ /pos/) || ($field =~/phys/)) {
		# Try to compare position fields by delta instead of absolute position
		if ($lastName eq "") {
		    # First marker, so go absolute
		    if ($map1{$name}{$field} != $map2{$name}{$field}) {
			$are_different += 1;
			if (($posFuzz == 0) || (abs($map1{$name}{$field} - $map2{$name}{$field}) > $posFuzz)) {
			    if ($tabular) {
				print "#A,$path1,$path2,$name,$field,".$map1{$name}{$field}.",".$map2{$name}{$field}.",".sprintf("%.5f", abs($map1{$name}{$field} - $map2{$name}{$field}))."\n";
			    } else {
				print "2: Marker \"$name\" has a different value for $field - 1:".$map1{$name}{$field}." vs 2:".$map2{$name}{$field}." (actual difference of ".
				    sprintf("%.5f", abs($map1{$name}{$field} - $map2{$name}{$field}))."\n";
			    }
			}
		    }
		} else {
		    # There is a previous marker, so use delta
		    if (($map1{$name}{$field} - $map1{$lastName}{$field}) != ($map2{$name}{$field} - $map2{$lastName}{$field})) {
			$are_different += 1;
			if (($posFuzz == 0) || (abs(($map1{$name}{$field} - $map1{$lastName}{$field}) - ($map2{$name}{$field} - $map2{$lastName}{$field})) > $posFuzz)) {
			    if ($tabular) {
				print "#R,$path1,$path2,$name,$lastName,$field,".$map1{$name}{$field}.",".$map1{$lastName}{$field}.",".$map2{$name}{$field}.",".$map2{$lastName}{$field}.",".
				    sprintf("%.5f", $map1{$name}{$field} - $map1{$lastName}{$field}).",".sprintf("%.5f", $map2{$name}{$field} - $map2{$lastName}{$field})."\n";
			    } else {
				print "2: Marker \"$name\" $field distance from predecessor \"$lastName\" is different - 1:".
				    sprintf("%.5f", $map1{$name}{$field} - $map1{$lastName}{$field})." vs 2: ".sprintf("%.5f", $map2{$name}{$field} - $map2{$lastName}{$field})."\n";
			    }
			}
		    }
		}
		$lastName = $name;
	    } else {
		if ($map1{$name}{$field} ne $map2{$name}{$field}) {
		    $are_different += 1;
		    print "2: Marker \"$name\" has a different value for $field - 1:".$map1{$name}{$field}." vs 2:".$map2{$name}{$field}."\n";
		}
	    }
	}
    }
}

if ($freqfiles) {
    if ($tabular) {
	print "#F,Path1,Path2,Marker,Allele,Freq1,Freq2,Difference\n";
    }
    # Describe and compare the allele frequency files
    # First build frequency file marker lists since they're not intrinsically present. A superset map file might be present.
    my @freqList1 = @{$$dataset1{maporder}};
    my %markers1 = %{$$dataset1{markers}};
    for (my $i = $#freqList1; $i >= 0; --$i) {
	splice(@freqList1, $i, 1) if (!defined($markers1{$freqList1[$i]}{alleles}));
    }
    my @type = ("", "microsatellite", "SNP", "microsatellite and SNP");
    print "1: Frequency file for ".($#freqList1 + 1)." ".$type[($$dataset1{snps}*2+$$dataset1{microsats})]." markers\n" if $verbose;
    my @freqList2 = @{$$dataset2{maporder}};
    my %markers2 = %{$$dataset2{markers}};
    for (my $i = $#freqList2; $i >= 0; --$i) {
	splice(@freqList2, $i, 1) if (!defined($markers2{$freqList2[$i]}{alleles}));
    }
    print "2: Frequency file for ".($#freqList2 + 1)." ".$type[($$dataset2{snps}*2+$$dataset2{microsats})]." markers\n" if $verbose;

    # Compare frequency files by looping over the superset marker names
    for my $name (uniqua (@freqList1, @freqList2)) {
	if (!defined($markers1{$name})) {
	    $are_different += 1;
	    print "1: Marker \"$name\" not found in frequency file, skipping!\n";
	    next;
	}
	if (!defined($markers2{$name})) {
	    $are_different += 1;
	    print "2: Marker \"$name\" not found in frequency file, skipping!\n";
	    next;
	}
	# Loop over superset of allele names
	my %alleles1 = %{$markers1{$name}{alleles}};
	my %alleles2 = %{$markers2{$name}{alleles}};
	for my $allele (uniqua (keys %alleles1, keys %alleles2)) {
	    if (!defined($alleles1{$allele})) {
		$are_different += 1;
		if ($tabular) {
		    print "#F,$path1,$path2,$name,$allele,N/A,".($alleles2{$allele}).",N/A\n";
		} else {
		    print "1: Marker \"$name\" allele \"$allele\" not found in frequency file, skipping!\n";
		}
		next;
	    }
	    if (!defined($alleles2{$allele})) {
		$are_different += 1;
		if ($tabular) {
		    print "#F,$path1,$path2,$name,$allele,".($alleles1{$allele}).",N/A,N/A\n";
		} else {
		    print "2: Marker \"$name\" allele \"$allele\" not found in frequency file, skipping!\n";
		}
		next;
	    }
	    if ($alleles1{$allele} != $alleles2{$allele}) {
		$are_different += 1;
		if (($freqFuzz == 0) || (abs($alleles1{$allele} - $alleles2{$allele}) > $freqFuzz)) {
		    if ($tabular) {
			print "#F,$path1,$path2,$name,$allele,".($alleles1{$allele}).",".($alleles2{$allele}).",".sprintf("%.5f", abs($alleles1{$allele} - $alleles2{$allele}))."\n";
		    } else {
			print "2: Marker \"$name\" allele \"$allele\" has a different frequency - 1:".
			    $alleles1{$allele}." vs 2:".$alleles2{$allele}." (actual difference of ".sprintf("%.5f", abs($alleles1{$allele} - $alleles2{$allele})).")\n";
		    }
		}
	    }
	}
    }
}

if ($pedfiles) {
    # Compare the superset of trait columns
    my %traits1 = %{$$dataset1{traits}};
    my %traits2 = %{$$dataset2{traits}};
    for my $trait (uniqua (keys %traits1, keys %traits2)) {
	if (!defined($traits1{$trait})) {
	    print "1: Trait \"$trait\" not found, skipping!\n";
	    $are_different += 1;
	    next;
	}
	if (!defined($traits2{$trait})) {
	    print "2: Trait \"$trait\" not found, skipping!\n";
	    $are_different += 1;
	    next;
	}
    }

    # First build locus lists since they're not intrinsically present. A superset map file might be present.
    my @locusList1 = @{$$dataset1{maporder}};
    my %markers1 = %{$$dataset1{markers}};
    for (my $i = $#locusList1; $i >= 0; --$i) {
	splice(@locusList1, $i, 1) if (!defined($markers1{$locusList1[$i]}{idx}));
    }
    my @locusList2 = @{$$dataset2{maporder}};
    my %markers2 = %{$$dataset2{markers}};
    for (my $i = $#locusList2; $i >= 0; --$i) {
	splice(@locusList2, $i, 1) if (!defined($markers2{$locusList2[$i]}{alleles}));
    }

    # Compare the superset of marker columns
    for my $marker (uniqua (@locusList1, @locusList2)) {
	if (!defined($markers1{$marker})) {
	    print "1: Marker \"$marker\" not found, skipping!\n";
	    $are_different += 1;
	    next;
	}
	if (!defined($markers2{$marker})) {
	    print "2: Marker \"$marker\" not found, skipping!\n";
	    $are_different += 1;
	    next;
	}
    }

    # Compare the pedigree files by looping over the superset of keys (pedids)
    for my $pedid (uniqua ((keys %families1, keys %families2))) {
	my %individuals1;
	my %individuals2;
	# Compare individuals by looping over the superset of keys (indids)
	if (!defined($families1{$pedid})) {
	    print "1: Pedigree \"$pedid\" not found, skipping!\n";
	    $are_different += 1;
	    next;
	}
	my @aref = @{$ {$families1{$pedid}}{individuals}};
	map { $individuals1{$_->indid} = $_; } @aref;

	if (!defined($families2{$pedid})) {
	    print "2: Pedigree \"$pedid\" not found, skipping!\n";
	    $are_different += 1;
	    next;
	}
	@aref = @{$ {$families2{$pedid}}{individuals}};
	map { $individuals2{$_->indid} = $_; } @aref;

	for my $indid (uniqua (keys %individuals1, keys %individuals2)) {
	    my %individual1;
	    my %individual2;
	    if (!defined($individuals1{$indid})) {
		print "1: Pedigree \"$pedid\" individual \"$indid\" not found, skipping!\n";
		$are_different += 1;
		next;
	    }
	    %individual1 = %{$individuals1{$indid}};
	    if (!defined($individuals2{$indid})) {
		print "2: Pedigree \"$pedid\" individual \"$indid\" not found, skipping!\n";
		$are_different += 1;
		next;
	    }
	    %individual2 = %{$individuals2{$indid}};

	    # Compare all scalar attributes (nice if they're named in an illuminating fashion)
	    for my $key (sort keys %individual1) {
		next if (($$dataset1{origfmt} eq "pre" or $$dataset2{origfmt} eq "pre") and
			 ($key eq "matsibid" or
			  $key eq "patsibid" or
			  $key eq "proband" or
			  $key eq "firstchildid"));
		next if ((!defined ($individual1{$key})) and (!defined ($individual2{$key})));
		next if (ref ($individual1{$key}) ne "");
		if ($individual1{$key} ne $individual2{$key}) {
		    $are_different += 1;
		    print "2: Ped \"$pedid\" ind \"$indid\" has different values for $key - 1:".$individual1{$key}." vs 2:".$individual2{$key}."\n";
		}
	    }

	    # Compare all common traits considering position
	    for my $trait (sort keys %traits1) {
		next if (!defined($traits2{$trait})); # We've already complained once, so just bail
		my $i = $traits1{$trait}{col}; my $j = $traits2{$trait}{col};
		if ($individual1{traits}[$i] ne $individual2{traits}[$j]) {
		    $are_different += 1;
		    print "2: Ped \"$pedid\" ind \"$indid\" has different values for trait \"$trait\" - 1:".
			$individual1{traits}[$i]." vs 2:".$individual2{traits}[$j]."\n";
		}
	    }

	    # Compare all common marker alleles for up to 10 markers per individual, otherwise swaps of SNP chip genotypes get crazy
	    my $markers_different = 0;
	    for my $locus (@locusList1) {
		next if (!defined($markers2{$locus}));
		my $i = $markers1{$locus}{idx}; my $j = $markers2{$locus}{idx};
		if (($individual1{markers}[$i][0] ne $individual2{markers}[$j][0]) or
		    ($individual1{markers}[$i][1] ne $individual2{markers}[$j][1])) {
		    $are_different += 1;
		    if ($markers_different++ < 10) {
			print "2: Ped \"$pedid\" (\"".$individual1{origpedid}."\") ind \"$indid\" (\"".$individual1{origindid}."\") has different values for marker ".($i+1)." (\"".$$dataset1{markerorder}[$i]."\") - 1:".
			    $individual1{markers}[$i][0]." ".$individual1{markers}[$i][1]." vs 2:".
			    $individual2{markers}[$j][0]." ".$individual2{markers}[$j][1]."\n";
		    }
		}
	    }
	    if ($markers_different > 0) {
		print "2: Ped \"$pedid\" (\"".$individual1{origpedid}."\") ind \"$indid\" (\"".$individual1{origindid}."\") has different values for $markers_different marker(s)\n";
	    }
	}
    }
}

warner ("Files are different") if ($are_different > 0);

sub uniqn {
    return sort { $a <=> $b } keys %{{ map { $_ => 1 } @_ }};
}

# Yeah, the Backyardigan!
sub uniqua {
    return sort keys %{{ map { $_ => 1 } @_ }};
}

#
# Check the paths to scripts we need to get work done. Lifted from Kelvin.
#
sub check_paths
{
    # For all of these, we allow environment variables to override everything,
    # even if values were set during installation.
    if ($ENV{KELVIN_ROOT}) {
        ($KELVIN_ROOT !~ /no_kelvin_root/i)
            and informer ("overriding installed KELVIN_ROOT with '$ENV{KELVIN_ROOT}' from environment");
        $KELVIN_ROOT = $ENV{KELVIN_ROOT};
    } elsif ($KELVIN_ROOT =~ /no_kelvin_root/i) {
        $KELVIN_ROOT = dirname ($0);
#        warner ("no KELVIN_ROOT defined by installation, using '$KELVIN_ROOT'");
    }

    # Instead of 'use'ing our support modules up above, we 'require' them now,
    # after we've had a chance to add KELVIN_ROOT to @INC.
    unshift (@INC, $KELVIN_ROOT);
    require KelvinDataset;
    require KelvinConfig;
    require KelvinFamily;
    # Check versions ('use' would have done for us automatically)
    KelvinDataset->VERSION (v1.3.0);
    KelvinConfig->VERSION (v1.2.0);
    KelvinFamily->VERSION (v1.3.0);
}

sub fatal
{   
    die ("FATAL - ABORTING, @_");
}

sub error
{   
    die ("ERROR - EXITING, @_\n");
}

sub warner
{   
    warn ("WARNING, @_\n");
}

sub informer
{   
    warn ("INFO, @_\n") if $verbose;
}

__END__


=head1 NAME

kdiff - validate and compare sets of Kelvin data files

=head1 SYNOPSIS

Use:

=over 5

kdiff [--verbose] [--pedfiles] [--freqfiles] [--mapfiles] [--c1 CONFIG1] [--c2 CONFIG2] [--p1 PATH1] [--p2 PATH2] [--freq-fuzz MAXDIFF] [--pos-fuzz MAXDIFF] [--tabular]

=back

where CONFIG1 and CONFIG2 are paths to standard Kelvin configuration files, and 
PATH1 and PATH2 are the default paths for locating data files referenced or defaulted
in CONFIG1 and CONFIG2.

=head1 DESCRIPTION

kdiff.pl validates and compares sets of Kelvin data files as identified by a pair of
(potentially dummy) configuration files. kdiff is primarily intended to validate the
transformations of data files performed as a part of the cleaning protocol. kdiff is
concerned only with content and not formatting or column order where it does not affect
the analysis.

If at least one of the data file options is specified, then only
the data files indicated by the options will be processed. If
no data file options are specified, then all data files
referenced (or defaulted) in the configuration file that actually exist
will be processed. As a reminder, the default Kelvin data files are "pedfile.dat",
"datafile.dat", "markers.dat" and "mapfile.dat".

Pedigree files require locus files for proper interpretation.

Configuration file(s) can be completely empty if the intention is to use the default 
Kelvin data file names, e.g.:

=over 5

kdiff --pedfiles --p1 old --p2 new --c1 /dev/null --c2 /dev/null

=back

will compare "old/pedfile.dat" using "old/datafile.dat" for column info
with "new/pedfile.dat" using "new/datafile.dat" for column info.

Config file(s) default to "./kelvin.conf", so:

=over 5

kdiff --p2 old

=back

will compare all data files referenced (or defaulted) in "./kelvin.conf" to
the identically-named data files in the "old" subdirectory. NOTE that "./kelvin.conf" is used
for both CONFIG1 and CONFIG2.  If you want to use the "old/kelvin.conf"
as CONFIG2, you need to specify it, e.g.:

=over 5

kdiff --p2 old --c2 old/kelvin.conf

=back

=head2 OPTIONS

=over 3

=item B<--verbose>

Provide extensive output describing the characteristics of
the data files in addition to information on differences.

=item B<--pedfiles>

Validate and compare pre- or post-makeped pedigree files as described 
by associated locus files. Validate locus files against 
frequency files if --freqfiles is specified as well. Validate 
locus files against map files if --mapfiles is specified as well.

=item B<--freqfiles>

Validate and compare frequency files. Validate frequency files
against locus files if --pedfiles is specified
as well. Validate frequency files against map files if --mapfile
is specified as well.

=item B<--mapfiles>

Validate and compare map files. Validate map files against locus
files if --pedfiles is specified as well. Validate
map files against frequency files if --freqfiles is specified as
well.

=item B<--c1 CONFIG1> | B<--configfile1 CONFIG1>

Specifies a standard Kelvin configuration file that references the data
files to be considered, or nothing at all if the default data file names
are to be used. If the option is not specified, CONFIG1 defaults to "./kelvin.conf".

If analysis characteristics are included in the configuration file, they will
be validated, and could cause kdiff to exit prematurely if incorrectly
specified.

=item B<--c1 CONFIG2> | B<--configfile1 CONFIG2>

A standard Kelvin configuration file. See B<--c1 CONFIG1> for details.

=item B<--p1 PATH1> | B<--path1 PATH1>

Locate data files referenced or defaulted in CONFIG1 relative to PATH1.
The current path is used if this option is not specified. Configuration
files are typically written with data file paths relative to
some presumed current default path for Kelvin execution. Since kdiff compares two
sets of data files, a single current default path cannot be used,
so B<--p1 PATH> and B<--p2 PATH> are provided to allow the specification of separate
default paths for each set of data files.

=item B<--p2 PATH2> | B<--path2 PATH2>

Locate data files referenced or defaulted in CONFIG1 relative to PATH2.
See B<--p1 PATH1> for details.

=item B<--freq-fuzz MAXDIFF>

Don't display any differing allele frequencies unless they are greater than MAXDIFF.
E.g. --freq 0.01 won't display frequencies of 0.1 and 0.999 as different, although 
the files will still be flagged as different.

=item B<--pos-fuzz MAXDIFF>

Don't display any different marker positions unless they are greater than MADIFF.
E.g. --pos 8 won't display positions of 87.234 and 93.500 as different, although
the files will still be flagged as different.

=item B<--tabular>

Display marker position and allele frequency differences in tabular (CSV) format. 
You'll want to grep them out of the rest of the output stream to use them. The
first column will be #F for allele frequencies, and #A and #R for marker positions.

Assuming you produce an output file called compare_all.txt from a mult-analysis run,
you can split-up the output with something like:

grep "^#A" compare_all.txt | sort -r | uniq >compare_all.absolute_positions.csv

grep "^#R" compare_all.txt | sort -r | uniq >compare_all.relative_positions.csv

grep "^#F" compare_all.txt | sort -r | uniq >compare_all.frequencies.csv

grep -v -e "^#A" -e "^#R" -e "^#F" compare_all.txt >compare_all.other.txt

=back

=head1 CAVEATS

This code is still very much under development.

=head1 COPYRIGHT

Copyright (C) 2011, 2022 Mathematical Medicine LLC
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.
You should have received a copy of the GNU General Public License along
with this program. If not, see <https://www.gnu.org/licenses/>.

=head1 AUTHOR

Bill Valentine-Cooper (William.Valentine-Cooper@NationwideChildrens.org)

=head1 DATE

$Id$

=cut
