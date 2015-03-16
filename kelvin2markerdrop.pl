#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min);
use Data::Dumper;
use File::Basename;
use lib dirname($0);
use kf;

#####################################
#
# $Id$
#
# by Bill Valentine-Cooper
#
# Copyright 2010, Nationwide Children's Research Institute.  All
# rights reserved.  Permission is hereby given to use this software
# for non-profit educational purposes only.
#

$| = 1;    # Force flush of output as printed.

# Sanctioned globals
my $PairCount = 0;
my $LiabilityClasses = 1;
my $maf0           = 0; # True if we found missing minor allele frequencies
my $pedFile          = "pedpost.dat";
my $mapFile          = "mapfile.dat";
my $companionFile    = "datafile.dat";
my $markersFile      = "markers.dat";                 # Default input files

# Heavy-duty diagnostics
my $DiagSolos = 0;    # Writes out a separate file for each marker with actual copies of template pedigrees

# Command line option flags
my $stats          = 0;
my $nokelvin       = 0;
my $write          = "unspecified";
my $WritePrefix    = "MD";

# Sorting routines
sub numerically { $a <=> $b }

sub numericIsh {
    if (($a . $b) =~ /^\d+$/) {
        $a <=> $b;
    } else {
        $a cmp $b;
    }
}

#####################################
#
sub dirPhenoCodes {
    $UnknownAffection = $Directives{PhenoCodes}[0];
    $Unaffected       = $Directives{PhenoCodes}[1];
    $Affected         = $Directives{PhenoCodes}[2];
}

#####################################
#
# Doesn't support much other than the type of simulation we currently care about.
#
# Verify command line parameters, check flags and do what the user asks.
#
print '$Id$'; print "\n";
my $Usage = <<EOF;

Usage "perl $0 [<flags>...] <kelvin configuration file>"

where <flags> are any of:

-nokelvin	Skip verification that kelvin can handle the analysis.
-stats		Print statistics on the make-up of the pedigree(s).
-write=[<prefix>] Write Morgan files, optionally with prefix <prefix> instead of "MD".

The configuration file will be read and analyzed, as well as the files that
it refers to.

Writes a Morgan format pedigree file and *template* parameter file that refers
to it. The parameter file is only a template because 

EOF

GetOptions(
    'nokelvin'       => \$nokelvin,
    'stats'          => \$stats,
    'write:s'        => \$write,
) or die "Invalid command line parameters.";
if ($write ne "unspecified") {
    $WritePrefix = $write if ($write ne "");
    $write = 1;
} else {
    $write = 0;
}

$Data::Dumper::Sortkeys = 1;

die "Invalid number of arguments supplied.\n$Usage"       if ($#ARGV < 0);
print "-nokelvin flag seen\n"                             if ($nokelvin);
print "-stats flag seen\n"                                if ($stats);
print "-write seen, using \"$WritePrefix\" prefix\n"      if ($write);

my $ConfFile = shift;
loadConf($ConfFile);
if (defined($Directives{LocusFile}[0])) { $companionFile = $Directives{LocusFile}[0]; }
loadCompanion($companionFile);
if (defined($Directives{MapFile}[0])) { $mapFile = $Directives{MapFile}[0]; }
loadMap($mapFile);
if (defined($Directives{FrequencyFile}[0])) { $markersFile = $Directives{FrequencyFile}[0]; }
loadFrequencies($markersFile);
if (defined($Directives{PedigreeFile}[0])) { $pedFile = $Directives{PedigreeFile}[0]; }

$maf0 = addMissingAlleles();

my $pedFileType = assessPedigree($pedFile, "POST", 0, 0);
$PairCount = loadPedigree($pedFile, $pedFileType, 0, 0);

#print Dumper(\%Pedigrees);
#print Dumper(\@Loci);
#print Dumper(\%LociAttributes);

checkRelations($pedFileType);
checkIntegrity();

if (!$nokelvin) {
    kelvinLimits();
}

if ($stats) {
    perfStats();
}

# Now for the meat of the script -- write a Morgan pedigree and markerdrop parameter file

open OUT, ">" . $WritePrefix . "_Pedigree.Dat";
print OUT '# $Id$'; print OUT "\n";

my $TotalInd = 0; # Morgan wants to know count of individuals before trying to read. Seems amateurish.
map { $TotalInd += scalar (keys (%{$Pedigrees{$_}})) } keys (%Pedigrees);

print OUT "input pedigree size ".$TotalInd."\n";
print OUT "input pedigree record names 3 integers 2\n\n";
print OUT "****************************************\n";

my $PedSeq = "a"; # This kind of craziness is what I love about Perl. I can increment this and go to multiple "digits" even.
# In case you're curious, I'm having to merge the individuals from all pedigrees into one giant disconnected "family" and 
# yet retain their individuality. We might be able to do this by prepending the pedigree number to individual numbers, but
# what if the combination overlaps with one from another pedigree, or the result is too wide? markerdrop only takes 8. This 
# handles both situations.
for my $Ped (sort numericIsh keys %Pedigrees) {
    for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
	my $Dad = $Pedigrees{$Ped}{$Ind}{Dad}; my $Mom = $Pedigrees{$Ped}{$Ind}{Mom};
	$Dad = $PedSeq.$Dad if ($Dad ne $UnknownPerson);
	$Mom = $PedSeq.$Mom if ($Mom ne $UnknownPerson);
	print OUT sprintf("%8s %8s %8s %1s %1s\n", $PedSeq.$Ind, $Dad, $Mom, $Pedigrees{$Ped}{$Ind}{Sex}, $Pedigrees{$Ped}{$Ind}{Aff});
    }
    $PedSeq++;
}
close OUT;

# Now for the markerdrop parameter file

open OUT, ">" . $WritePrefix . "_Parameter.Dat";
print OUT '# $Id$'; print OUT "\n";

print OUT "simulate markers ".(scalar(@Loci))." using trait\n";
print OUT "map marker Kosambi positions ";

for my $Locus (@Loci) {
    next if ($LociAttributes{$Locus}{Type} eq "T");
    print OUT $Map{$Locus}{SAPos}." ";
}
print OUT "\n";

print OUT "map trait 1 marker <nth> dist <cM>\n";
print OUT "set incomplete penetrance <dd> <dD and Dd> <DD>\n";
print OUT "set trait 1 freqs <d> <D>\n";
print OUT "input seed file './trait_and_marker.seed'\n";
print OUT "output overwrite seed file './trait_and_marker.seed'\n";
print OUT "input pedigree file './".$WritePrefix."_Pedigree.Dat'\n";
print OUT "output overwrite pedigree file './".$WritePrefix."_Pedigree.Out'\n";

my $MarkerSeq = 0;
for my $Locus (@Loci) {
    next if ($LociAttributes{$Locus}{Type} eq "T");
    print OUT "set markers ".++$MarkerSeq." freqs ";
    for my $Allele (@{$LociAttributes{ $Locus }{Alleles}{OrderedList}}) {
        print "Marker $Locus allele $Allele\n"; # has frequency ".$LociAttributes{ $Locus }{Alleles}{$Allele}{Frequency}."\n";
        print OUT $LociAttributes{ $Locus }{Alleles}{$Allele}{Frequency}." ";
    }
    print OUT "\n";
}

print OUT "\nset markers $MarkerSeq data\n";

$PedSeq = "a";
for my $Ped (sort numericIsh keys %Pedigrees) {
    for my $Ind (sort numericIsh keys %{ $Pedigrees{$Ped} }) {
        print OUT $PedSeq.$Ind;
        if ( shift ( @{ $Pedigrees{$Ped}{$Ind}{Mks} } ) ne "0 0") {
            print OUT ("  1 1") x $MarkerSeq;
        } else {
            print OUT ("  0 0") x $MarkerSeq;
        }
        print OUT "\n";
    }
    $PedSeq++;
}

close OUT;

exit;
