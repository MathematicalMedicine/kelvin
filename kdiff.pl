#!/usr/bin/perl -w

use strict;
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
my $freqfiles = 0; my $dup_freqfiles = 0;
my $mapfiles = 0; my $dup_mapfiles = 0;
my $path1 = ''; my $path2 = '';
my $pedfile1 = ''; my $pedfile2 = '';
my $locusfile1 = ''; my $locusfile2 = '';
my $freqfile1 = ''; my $freqfile2 = '';
my $mapfile1 = ''; my $mapfile2 = '';
my $configfile1 = './kelvin.conf'; my $configfile2 = './kelvin.conf';

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
    'path2|p2=s' => \$path2
    ) or pod2usage(2);

pod2usage(-noperldoc => 1, -verbose => 2) if ($help); # All source is shown if -noperdoc isn't specified!

pod2usage(-verbose => 0) if ($#ARGV > 1);

check_paths ();

print "Processing configuration files 1:\"$configfile1\" and 2:\"$configfile2\"\n" if $verbose;

check_paths ();

$use_found_files = 1 if (($pedfiles eq 0) and ($freqfiles eq 0) and ($mapfiles eq 0));
print "Processing all \"found\" data files\n" if ($verbose and $use_found_files);

# First we load (and thereby validate) everything

($config1 = KelvinConfig->new ($configfile1))
    or error ("Processing of configuration file 1:\"$configfile1\" failed: $KelvinConfig::errstr");
($config2 = KelvinConfig->new ($configfile2))
    or error ("Processing of configuration file 2:\"$configfile2\" failed: $KelvinConfig::errstr");

if ($use_found_files or $pedfiles) {
    # Construct path and name for pedigree and locus files
    $$data1ref{PedigreeFile} = $ {$config1->isConfigured ("PedigreeFile")}[0];
    $$data1ref{LocusFile} = $ {$config1->isConfigured ("LocusFile")}[0];
    # Consider any specified path1..
    $$data1ref{PedigreeFile} = catfile($path1, $$data1ref{PedigreeFile})
	if ((!file_name_is_absolute($$data1ref{PedigreeFile})) && $path1);
    $$data1ref{LocusFile} = catfile($path1, $$data1ref{LocusFile})
	if ((!file_name_is_absolute($$data1ref{LocusFile})) && $path1);

    $$data2ref{PedigreeFile} = $ {$config2->isConfigured ("PedigreeFile")}[0];
    $$data2ref{LocusFile} = $ {$config2->isConfigured ("LocusFile")}[0];
    # Consider any specified path2...
    $$data2ref{PedigreeFile} = catfile($path2, $$data2ref{PedigreeFile})
	if ((!file_name_is_absolute($$data2ref{PedigreeFile})) && $path2);
    $$data2ref{LocusFile} = catfile($path2, $$data2ref{LocusFile})
	if ((!file_name_is_absolute($$data2ref{LocusFile})) && $path2);

    # Turn on option if we find them, but don't turn it off if we don't so we can error-out as needed.
    $pedfiles = 1 if ((-e $$data1ref{PedigreeFile}) and (-e $$data1ref{LocusFile}) and
		      (-e $$data2ref{PedigreeFile}) and (-e $$data2ref{LocusFile}));

    # If all files are the same, don't bother validating again. We'll still compare just to verify this code.
    $dup_pedfiles = 1 if ((($$data1ref{PedigreeFile}) eq ($$data2ref{PedigreeFile})) and
			  (($$data1ref{LocusFile}) eq ($$data2ref{LocusFile})));
}

if ($pedfiles) {
    # Read and validate the pedigree and locus files
    print "Validating pedigree file 1:\"".$$data1ref{PedigreeFile}."\" as described by locus file 1:\"".$$data1ref{LocusFile}."\"\n" if $verbose;
    $dataset1 = KelvinDataset->new ($data1ref)
	or error ("Processing of pedigree file 1:\"".$$data1ref{PedigreeFile}."\" as described by locus file 1:\"".$$data1ref{LocusFile}."\ failed: $KelvinDataset::errstr");
    my $family;
    my $totalInds = 0;
    # print "Dataset structure is ".Dumper($dataset1)."\n";
    while ($family = $dataset1->readFamily) {
	# print "Family structure is ".Dumper($family)."\n";
	$$dataset1{origfmt} = $$family{origfmt};
	$totalInds += $$family{count};
	print $family->pedtype." family ".$$family{pedid}." of ".$$family{count}." (".$$family{founders}."f/".$$family{nonfounders}."nf)\n" if $verbose;
	$families1{$$family{pedid}} = $family;
    }
    (defined ($family))
	or error ("Read of pedigree in file 1 failed, $KelvinDataset::errstr");
    print "P/L 1: A ".$$dataset1{origfmt}."-makeped format file with ".scalar(keys %{$$dataset1{markers}})." markers, ".
	scalar(keys %families1)." families and $totalInds individuals.\n";

    if ($dup_pedfiles) {
	$dataset2 = $dataset1;
	%families2 = %families1;
    } else {
	print "Validating pedigree file 2:\"".$$data2ref{PedigreeFile}."\" as described by locus file 2:\"".$$data2ref{LocusFile}."\"\n" if $verbose;
	$dataset2 = KelvinDataset->new ($data2ref)
	    or error ("Processing of pedigree file 2:\"".$$data2ref{PedigreeFile}."\" as described by locus file 2:\"".$$data2ref{LocusFile}."\ failed: $KelvinDataset::errstr");
	$totalInds = 0;
	while ($family = $dataset2->readFamily) {
	    # print "Family structure is ".Dumper($family)."\n";
	    $$dataset2{origfmt} = $$family{origfmt};
	    $totalInds += $$family{count};
	    print $family->pedtype." family ".$$family{pedid}." of ".$$family{count}." (".$$family{founders}."f/".$$family{nonfounders}."nf)\n" if $verbose;
	    $families2{$$family{pedid}} = $family;
	}
	(defined ($family))
	    or error ("Read of pedigree in file 2 failed, $KelvinDataset::errstr");
	print "P/L 2: A ".$$dataset2{origfmt}."-makeped format file with ".scalar(keys %{$$dataset2{markers}})." markers, ".
	    scalar(keys %families2)." families and $totalInds individuals.\n";
    }
}

if ($freqfiles) {
    # Read and validate the frequency files
    print "WARNING -- Read/validate frequency files not implemented yet.";
    exit 1;
}

if ($mapfiles) {
    # Read and validate the map files
    print "WARNING -- Read/validate map files not implemented yet.";
    exit 1;
}

if ($pedfiles) {
    # Compare the pedigree files by looping over the superset of keys (pedids)
    for my $pedid (uniq ((keys %families1, keys %families2))) {
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

	for my $indid (uniq (keys %individuals1, keys %individuals2)) {
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
	    for my $key (keys %individual1) {
		next if ((!defined ($individual1{$key})) and (!defined ($individual2{$key})));
		next if (ref ($individual1{$key}) ne "");
		if ($individual1{$key} ne $individual2{$key}) {
		    $are_different += 1;
		    print "Ped \"$pedid\" ind \"$indid\" has different values for $key - 1:".$individual1{$key}." vs 2:".$individual2{$key}."\n";
		}
	    }

	    # All traits...(every one!)
	    for (my $i=0; $i<1; $i++) {
		if ($individual1{traits}[$i] ne $individual2{traits}[$i]) {
		    $are_different += 1;
		    print "Ped \"$pedid\" ind \"$indid\" has different values for trait ".($i+1)." - 1:".
			$individual1{traits}[$i]." vs 2:".$individual2{traits}[$i]."\n";
		}
	    }

	    # All markers...
	    for (my $i=0; $i<scalar(keys %{$$dataset1{markers}}); $i++) {
		if (($individual1{markers}[$i][0] ne $individual2{markers}[$i][0]) or
		    ($individual1{markers}[$i][1] ne $individual2{markers}[$i][1])) {
		    $are_different += 1;
		    print "Ped \"$pedid\" ind \"$indid\" has different values for marker ".($i+1)." (\"".$$dataset1{markerorder}[$i]."\") - 1:".
			$individual1{markers}[$i][0]." ".$individual1{markers}[$i][1]." vs 2:".
			$individual2{markers}[$i][0]." ".$individual2{markers}[$i][1]."\n";
		}
	    }
	}
    }
    error ("Files are different") if ($are_different > 0);
}

sub uniq {
    return sort { $a <=> $b } keys %{{ map { $_ => 1 } @_ }};
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
            and warner ("overriding installed KELVIN_ROOT with '$ENV{KELVIN_ROOT}' from environment");
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
    KelvinDataset->VERSION (1.30);
    KelvinConfig->VERSION (1.20);
    KelvinFamily->VERSION (1.30);
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

__END__


=head1 NAME

kdiff - validate and compare sets of Kelvin data files

=head1 SYNOPSIS

Use:

=over 5

kdiff [--verbose] [--pedfiles] [--freqfiles] [--mapfiles] [--c1 CONFIG1] [--c2 CONFIG2] [--p1 PATH1] [--p2 PATH2] 

where CONFIG1 and CONFIG2 are paths to standard Kelvin configuration files, and 
PATH1 and PATH2 are the default paths for locating data files referenced or defaulted
in CONFIG1 and CONFIG2.

=back

=head1 DESCRIPTION

kdiff.pl validates and compares sets of Kelvin data files as identified by a pair of
(potentially dummy) configuration files. kdiff is primarily intended to validate the
transformations of data files performed as a part of the cleaning protocol.

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

=back

=head1 CAVEATS

This code is still very much under development.

=head1 COPYRIGHT

Copyright 2011, Nationwide Children's Hospital Research Institute
All rights reserved. Permission is hereby granted to use this software
for non-profit educational purposes only.

=head1 DATE

$Id$

=cut
