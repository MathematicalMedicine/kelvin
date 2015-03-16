#!/usr/bin/env perl -w

use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions;
use Data::Dumper;
no warnings 'once';

$|=1; # Immediate output

# kinfo --help for documentation.

# Used globally
my $svn_version='$Id$';

# Argument- and option-related
my $verbose = 0;
my $help = 0;
my $pedigrees = 1;
my $markers = 1;
my $map = 1;
my $path = '';
my $pedfile = '';
my $locusfile = '';
my $freqfile = '';
my $mapfile = '';
my $configfile = './kelvin.conf';

# Hash references to keep KelvinDataset happy
my ($dataref);

# Primary Kelvin* workhorses
my ($config);
my ($dataset);

# Something I'd expected to find in Kelvin* but wasn't there...
my (%families);

my $KELVIN_ROOT='no_kelvin_root';

GetOptions (
    'verbose' => \$verbose,
    'help|?' => \$help,
    'pedigrees' => \$pedigrees,
    'markers' => \$markers,
    'map' => \$map,
    'configfile|c=s' => \$configfile,
    'path|p=s' => \$path,
    ) or pod2usage(2);

pod2usage(-noperldoc => 1, -verbose => 2) if ($help); # All source is shown if -noperdoc isn't specified!

pod2usage(-verbose => 0) if ($#ARGV > -1);

check_paths ();

print "Processing configuration file:\"$configfile\"\n";

check_paths ();

# First we load (and thereby validate) everything

($config = KelvinConfig->new ($configfile))
    or error ("Processing of configuration file \"$configfile\" failed: $KelvinConfig::errstr");
my %phenocodes;
my %phenomap;
my ($UNDEFINED, $UNAFFECTED, $AFFECTED) = (0, 1, 2);
my $aref = $config->isConfigured ("PhenoCodes");
@phenocodes{$UNDEFINED, $UNAFFECTED, $AFFECTED} = split (/,\s*/, $$aref[0]);
map { defined ($phenocodes{$_}) and $phenomap{$phenocodes{$_}} = $_ } keys (%phenocodes);

if ($markers) {
    # Construct path and name for the frequency file using default since KelvinConfig doesn't.
    my $filename = defined ($config->isConfigured ("FrequencyFile")) ?
	$ {$config->isConfigured ("FrequencyFile")}[0] : "markers.dat";

    $filename = catfile($path, $filename)
	if ((!file_name_is_absolute($filename)) && $path);
    # Admitted strangeness here...I'm letting KelvinConfig give them the bad news.
    $$dataref{FrequencyFile} = $filename if ((-e $filename) or $markers);
}

if ($map) {
    # Construct path and name for the map file
    my $filename = $ {$config->isConfigured ("MapFile")}[0];
    $filename = catfile($path, $filename)
	if ((!file_name_is_absolute($filename)) && $path);
    $$dataref{MapFile} = $filename if ((-e $filename) or $map);
}

if ($pedigrees) {
    # Construct path and name for pedigree and locus files
    my $pedname = $ {$config->isConfigured ("PedigreeFile")}[0];
    my $locusname = $ {$config->isConfigured ("LocusFile")}[0];
    $pedname = catfile($path, $pedname)
	if ((!file_name_is_absolute($pedname)) && $path);
    $locusname = catfile($path, $locusname)
	if ((!file_name_is_absolute($locusname)) && $path);
    if (((-e $pedname) and (-e $locusname)) or $pedigrees) {
	$$dataref{PedigreeFile} = $pedname;
	$$dataref{LocusFile} = $locusname;
    }
}

# Read and validate everything we have (except the pedfile contents)
print "Validating files: ";
for my $key (sort %{$dataref}) {print "$key->".$$dataref{$key}." " if (defined($$dataref{$key})); }
print "\n";

$dataset = KelvinDataset->new ($dataref)
    or error ("Validation of referenced or defaulted data files failed: $KelvinDataset::errstr");
$dataset->setUndefPhenocode ($UNDEFINED);
my $markerOrder = $dataset->markerOrder;

if ($pedigrees) {
    # Read, validate and describe the pedigrees
    my $family; my $totalInds;

    my $pc=5; # Five percent
    my $mc=(scalar(@{$$dataset{maporder}}) + 1) * $pc / 100;

    print "Ped\tType\t#Ind\t#Fdr\t#Nfdr\t2N-F\t#Gen>$pc\%\t#Phen\t#Aff\t#G+A\t#Dum\t%Het\t%Msg\n";
    while ($family = $dataset->readFamily) {
	$$dataset{origfmt} = $$family{origfmt};
	$totalInds += $$family{count};
	my ($gC, $pC, $aC, $gaC, $dmC) = (0, 0, 0, 0, 0);
	my @g = (); # List of individuals with a genotype for finding degenerate cases
	my $het = 0; my $hom = 0; my $msgmk = 0; my $totmk = 0;
	foreach my $ind (@{$family->individuals}) {
	    if ($ind->{genotyped} > $mc) {
		$gC++;
		push @g, $ind;
	    }
	    $aC++ if (defined($ind->{traits}[0]) and ($ind->{traits}[0] == $AFFECTED));
	    $pC++ if (defined($ind->{phenotyped}) and ($ind->{phenotyped} > 0));
	    $gaC++ if ($ind->{genotyped} > $mc and (defined($ind->{traits}[0]) and ($ind->{traits}[0] == $AFFECTED)));
	    $dmC++ if ($ind->{genotyped} <= $mc and ($ind->{phenotyped} == 0));

	    # Determine heterozygosity and missingness
	    foreach my $marker (@$markerOrder) {
		my $aref = $ind->getGenotype ($marker);
		$totmk++;
		if ($$aref[0] eq "0") {
		    $msgmk++;
		    next;
		}
		# print "$$aref[0]/$$aref[1], ";
		($$aref[0] ne $$aref[1]) and $het++;
		($$aref[0] eq $$aref[1]) and $hom++;
	    }
	}
	# Identify degenerate cases IN A NON-DESTRUCTIVE MANNER for later removal
	# If the genotyped individuals are a parent and child, call this a uninformative duo. If these are two parents and 
	# a single child, call this an uninformative trio. Believe it or nuts this performs the fewest tests
	# to determine both cases.
	my $f=$family->pedtype;
	if ($gC == 2 and
	    ($g[0]->{momid} eq $g[1]->{indid} or $g[0]->{dadid} eq $g[1]->{indid} or
	     $g[1]->{momid} eq $g[0]->{indid} or $g[1]->{dadid} eq $g[0]->{indid})) {
	    $f = "ui-duo";
	} elsif ($gC == 3 and
		 ((($g[0]->{momid} eq $g[1]->{indid} and $g[0]->{dadid} eq $g[2]->{indid}) or
		   ($g[0]->{momid} eq $g[2]->{indid} and $g[0]->{dadid} eq $g[1]->{indid})) or

		  (($g[1]->{momid} eq $g[0]->{indid} and $g[1]->{dadid} eq $g[2]->{indid}) or
		   ($g[1]->{momid} eq $g[2]->{indid} and $g[1]->{dadid} eq $g[0]->{indid})) or

		  (($g[2]->{momid} eq $g[1]->{indid} and $g[2]->{dadid} eq $g[0]->{indid}) or
		   ($g[2]->{momid} eq $g[0]->{indid} and $g[2]->{dadid} eq $g[1]->{indid})))) {
	    $f = "ui-trio";
	} elsif ($gC == 1) {
	    $f = "ui-sglt";
	} elsif ($gC == 0) {
	    $f = "ui-null";
	}
	printf $$family{pedid}."\t".sprintf("%7.7s\t%d\t%d\t%d\t%d\t",$f, $$family{count}, $$family{founders}, $$family{nonfounders}, (2 * $$family{nonfounders}) - $$family{founders});
	print $gC."\t".$pC."\t".$aC."\t".$gaC."\t".$dmC."\t".int($het*100/($totmk > 0 ? $totmk : 1))."\t".int($msgmk*100/($totmk > 0 ? $totmk : 1))."\n";
	$families{$$family{pedid}} = $family;
    }
    (defined ($family))
	or error ("Read of pedigree file \"".$$dataref{PedigreeFile}."\" failed, $KelvinDataset::errstr");
    print $$dataset{origfmt}."-makeped format file with ".scalar(@{$$dataset{markerorder}})." markers, ".
	scalar(keys %families)." families and $totalInds individuals.\n";
}

if ($map) {
    # Describe the map file
    print $$dataset{mapfunction}." map of ".scalar(@{$$dataset{maporder}})." markers for chromosome ".$$dataset{chromosome}." providing ".join(", ",@{$$dataset{mapfields}})."\n";

    my %fieldhash = map { $_ => 1 } @{$$dataset{mapfields}};
}

if ($markers) {
    # Describe the allele frequency file
    # First build frequency file marker list since they're not intrinsically present. A superset map file might be present.
    my @freqList = @{$$dataset{maporder}};
    my %markers = %{$$dataset{markers}};
    for (my $i = $#freqList; $i >= 0; --$i) {
	splice(@freqList, $i, 1) if (!defined($markers{$freqList[$i]}{alleles}));
    }
    my @type = ("", "microsatellite", "SNP", "microsatellite and SNP");
    print "Frequency file for ".($#freqList + 1)." ".$type[($$dataset{snps}*2+$$dataset{microsats})]." markers\n";
}

if ($pedigrees) {
    # Describe trait column(s)
    my %traits = %{$$dataset{traits}};

    # First build locus lists since they're not intrinsically present. A superset map file might be present.
    my @locusList = @{$$dataset{maporder}};
    my %markers = %{$$dataset{markers}};
    for (my $i = $#locusList; $i >= 0; --$i) {
	splice(@locusList, $i, 1) if (!defined($markers{$locusList[$i]}{idx}));
    }
}

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

sub informer
{   
    warn ("INFO, @_\n") if $verbose;
}

__END__


=head1 NAME

kinfo - validate and describe Kelvin data files

=head1 SYNOPSIS

Use:

=over 5

kinfo [--verbose] [--pedigrees] [--markers] [--map] [--c CONFIG] [--p PATH]

=back

where CONFIG is a standard Kelvin configuration file, and PATH is the default path
 for locating data files referenced or defaulted in CONFIG.

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


=head1 CAVEATS

This code is still very much under development.

=head1 COPYRIGHT

Copyright 2011, Nationwide Children's Hospital Research Institute
All rights reserved. Permission is hereby granted to use this software
for non-profit educational purposes only.

=head1 AUTHOR

Bill Valentine-Cooper (William.Valentine-Cooper@NationwideChildrens.org)

=head1 DATE

$Id$

=cut
