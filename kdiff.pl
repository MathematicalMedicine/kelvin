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

my $verbose = 0;
my $help = 0;
my $pedfiles = 0;
my $freqfiles = 0;
my $mapfiles = 0;
my $path1 = '';
my $path2 = '';

my $pedfile1 = '';
my $pedfile2 = '';
my $locusfile1 = '';
my $locusfile2 = '';
my $freqfile1 = '';
my $freqfile2 = '';
my $mapfile1 = '';
my $mapfile2 = '';

my $config1;
my $config2;

# Hash references to keep KelvinDataset happy
my $data1ref;
my $data2ref;

my $dataset1;
my $dataset2;

my %families1;
my %families2;

my $KELVIN_ROOT='no_kelvin_root';

GetOptions (
    'verbose' => \$verbose,
    'help|?' => \$help,
    'pedfiles' => \$pedfiles,
    'freqfiles' => \$freqfiles,
    'mapfiles' => \$mapfiles,
    'path1=s' => \$path1,
    'path2=s' => \$path2
    ) or pod2usage(2);
pod2usage(1) if ($help);

pod2usage(2) if ($#ARGV != 1);

check_paths ();

my $configfile1 = shift();
my $configfile2 = shift();

print "Working with configuration files 1:\"$configfile1\" and 2:\"$configfile2\"\n" if $verbose;

check_paths ();

# First we load (and thereby validate) everything

($config1 = KelvinConfig->new ($configfile1))
    or error ("Processing of configuration file 1:\"$configfile1\" failed: $KelvinConfig::errstr");
($config2 = KelvinConfig->new ($configfile2))
    or error ("Processing of configuration file 2:\"$configfile2\" failed: $KelvinConfig::errstr");

if ($pedfiles) {
    # Read and validate the pedigree and locus files
    $$data1ref{PedigreeFile} = $ {$config1->isConfigured ("PedigreeFile")}[0];
    $$data1ref{LocusFile} = ${$config1->isConfigured ("LocusFile")}[0];
    # Consider any specified path1..
    $$data1ref{PedigreeFile} = catfile($path1, $$data1ref{PedigreeFile})
        if ((!file_name_is_absolute($$data1ref{PedigreeFile})) && $path1);
    $$data1ref{LocusFile} = catfile($path1, $$data1ref{LocusFile})
	if ((!file_name_is_absolute($$data1ref{LocusFile})) && $path1);
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
    print "P/D 1: A ".$$dataset1{origfmt}."-makeped format file with ".scalar(keys %{$$dataset1{markers}})." markers, ".
	scalar(keys %families1)." families and $totalInds individuals.\n";

    $$data2ref{PedigreeFile} = $ {$config2->isConfigured ("PedigreeFile")}[0];
    $$data2ref{LocusFile} = $ {$config2->isConfigured ("LocusFile")}[0];
    # Consider any specified path2...
    $$data2ref{PedigreeFile} = catfile($path2, $$data2ref{PedigreeFile})
	if ((!file_name_is_absolute($$data2ref{PedigreeFile})) && $path2);
    $$data2ref{LocusFile} = catfile($path2, $$data2ref{LocusFile})
	if ((!file_name_is_absolute($$data2ref{LocusFile})) && $path2);
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
    print "P/D 2: A ".$$dataset2{origfmt}."-makeped format file with ".scalar(keys %{$$dataset2{markers}})." markers, ".
	scalar(keys %families2)." families and $totalInds individuals.\n";

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
	    print "1: Pedigree $pedid not found, skipping!\n";
	    next;
	}
	my @aref = @{${$families1{$pedid}}{individuals}};
	map { $individuals1{$_->indid} = $_; } @aref;
	if (!defined($families2{$pedid})) {
	    print "2: Pedigree $pedid not found, skipping!\n";
	    next;
	}
	@aref = @{${$families2{$pedid}}{individuals}};
	map { $individuals2{$_->indid} = $_; } @aref;

	for my $indid (uniq (keys %individuals1, keys %individuals2)) {
	    my %individual1;
	    my %individual2;
	    if (!defined($individuals1{$indid})) {
		print "1: Pedigree $pedid individual $indid not found, skipping!\n";
		next;
	    }
	    %individual1 = %{$individuals1{$indid}};
	    if (!defined($individuals2{$indid})) {
		print "2: Pedigree $pedid individual $indid not found, skipping!\n";
		next;
	    }
	    %individual2 = %{$individuals2{$indid}};

	    # Compare all scalar attributes (nice if they're named in an illuminating fashion)
	    for my $key (keys %individual1) {
		next if ((!defined ($individual1{$key})) and (!defined ($individual2{$key})));
		next if (ref ($individual1{$key}) ne "");
		print "Pedigree $pedid individual $indid has different values for $key - 1:".$individual1{$key}." vs 2:".$individual2{$key}."\n"
		    if ($individual1{$key} ne $individual2{$key});
	    }

	    # All traits...(every one!)
	    for (my $i=0; $i<1; $i++) {
		print "Pedigree $pedid individual $indid has different values for trait $i - 1:".
		    $individual1{traits}[$i]." vs 2:".$individual2{traits}[$i]."\n" if ($individual1{traits}[$i] ne $individual2{traits}[$i]);
	    }

	    # All markers...
	    for (my $i=0; $i<scalar(keys %{$$dataset1{markers}}); $i++) {
		print "Pedigree $pedid individual $indid has different values for marker ".($i+1)." (".$$dataset1{markerorder}[$i].") - 1:".
		    $individual1{markers}[$i][0]." ".$individual1{markers}[$i][1]." vs 2:".
		    $individual2{markers}[$i][0]." ".$individual2{markers}[$i][1]."\n"
		    if (($individual1{markers}[$i][0] ne $individual2{markers}[$i][0]) or
			($individual1{markers}[$i][1] ne $individual2{markers}[$i][1]));
	    }
	}
    }
    
}

sub uniq {
    return sort { $a <=> $b } keys %{{ map { $_ => 1 } @_ }};
}

#
# Check the paths to scripts we need to get work done.
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

kdiff - validate and compare sets of kelvin data files

=head1 SYNOPSIS

kdiff [--verbose] [--pedfiles] [--freqfiles] [--mapfiles] [--path1] [--path2] configfile1 configfile2

=head1 DESCRIPTION

kdiff.pl validates and compares sets of kelvin data files as identified by a pair of
(potentially dummy) configuration files. kdiff is primarily intended to validate the
transformations of data files performed by cleaning.

At least one of the data file comparison flags needs to be specified.
Pedigree files require locus files for proper interpretation. Config
file(s) can be empty if Kelvin default data file names are used and
path is specified, e.g.:

kdiff --postfiles --path1 old --path2 new /dev/null /dev/null

will compare old/pedfile.dat using old/datafile.dat for column info
with new/pedfile.dat using new/datafile.dat for column info.

=head2 ARGUMENTS

=item B<configfile1 configfile1>

Standard Kelvin configuration files that need only specify the data
files to be considered, or nothing at all if the default data file names
are used and a path option is specified. Note that if analysis
characteristics are specified in the configuration files, they will
be validated, and could therefore cause kdiff to exit if incorrectly
specified.

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

=item B<--path1>

Locate data files specified in B<configfile1> relative to B<path1>.
Current path is used if B<--path1> is not specified. Configuration
files are typically written with data file paths relative to
some presumed current default path. Since kdiff compares two
sets of data files, a single current default path cannot be used,
so B<--path1> and B<--path2> are provided to allow the specification of a
"default path" for each set of data files instead.

=item B<--path2>

Locate data files specified in B<configfile2> relative to B<path2>.
See B<--path1> for details.

=back

=head1 CAVEATS

This code is still very much under development.

=head1 COPYRIGHT

Copyright 2011, Nationwide Children's Hospital Research Institute
All rights reserved. Permission is hereby granted to use this software
for non-profit educational purposes only.

=cut

