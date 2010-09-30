#!perl -w
use strict;
use KelvinConfig;
use constant {
    INDICATOR => 0,
    NAME => 1,
    CHR => 0,
    POS => 1,
    PPL => 2,
    BR => 3,
    LIST => 4
    };

# read locus file, create list of markers
# read BR file, create list of positions and associated markers

# create merlin directory
# create pre-makeped ped file
# copy freq and map files
# create one locus file for each position in the BR file

my $config;
my ($pedfile, $locusfile, $mapfile, $freqfile, $brfile);

my ($line, $lineno);
my $locilist = '';
my @fields;
my @groups = ();
my @positions = ();
my @traits = ();
my @markers = ();
my @loci;


if ($ARGV[0]) {
    $config = KelvinConfig->new ($ARGV[0]);
} elsif (-d "kelvin.conf") {
    $config = KelvinConfig->new ("kelvin.conf");
} else {
    die ("usage: $0 <kelvinconfigfile>\n");
}
($config) or die ("new KelvinConfig failed, $KelvinConfig::errstr\n");
$pedfile = $ {$config->isConfigured ("PedigreeFile")}[0];
$locusfile = $ {$config->isConfigured ("LocusFile")}[0];
$mapfile = $ {$config->isConfigured ("MapFile")}[0];
$freqfile = $ {$config->isConfigured ("FrequencyFile")}[0];
$brfile = $ {$config->isConfigured ("BayesRatioFile")}[0] . "-baseline";
(-f $brfile) or die ("expected a baseline BR file named '$brfile', but didn't find one\n");

open (IN, $brfile) or die ("open '$brfile failed, $!\n");
$lineno = 0;
while ($line = <IN>) {
    $lineno++;
    $line =~ /(\#|Chr)/ and next;
    @fields = ($line =~ /\s*(\d+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.e\+\-]+)\s+\(([\d\,]+)\)/)
	or die ("$brfile, line $lineno is badly formed\n");
    $fields[LIST] =~ s/,/+/;
    if ($locilist ne $fields[LIST]) {
	$locilist = $fields[LIST];
	push (@groups, $locilist);
	push (@positions, [ $fields[POS] ]);
    } else {
	push (@{$positions[-1]}, $fields[POS]);
    }
}

open (IN, $locusfile) or die ("open '$locusfile' failed, $!\n");
$lineno = 0;
while ($line = <IN>) {
    $lineno++;
    @fields = ($line =~ /\s*([ATCM])\s+(\w+)/)
	or die ("$locusfile, line $lineno is badly formed\n");
    if ($fields[INDICATOR] eq 'M') {
	push (@markers, $fields[NAME]);
    } else { 
	# $indicator is T, A or C
	$fields[INDICATOR] =~ s/T/A/;
	push (@traits, "$fields[INDICATOR] $fields[NAME]");
    }
}
close (IN);

print (join (',', @groups), "\n");

foreach $locilist (@groups) {
    @loci = split (/\+/, $locilist);
    foreach (@loci) { $_--; }

    if (! -d "markers_$locilist") {
	mkdir ("markers_$locilist") or die ("mkdir 'markers_$locilist' failed, $!\n");
    }
    open (OUT, ">markers_$locilist/positions")
	or die ("open 'markers_$locilist/positions' failed, $!\n");
    print (OUT map { "$_\n" } @{shift (@positions)});
    system ("qtdt_munge.pl", "-d", $locusfile, "-f", $freqfile, "-m", $mapfile, "-p", $pedfile,
	    "--subset",  "$markers[$loci[0]],$markers[$loci[-1]]", "--postmakeped", "--merlin",
	    "--prefix", "markers_$locilist/merlin");
    ($? != 0) and die ("qtdt_munge.pl failed\n");
    
    # Convert pedfile to pre-MAKEPED 
    open (IN, "markers_$locilist/merlin.post")
	or die ("open 'markers_$locilist/merlin.post' failed, $!\n");
    open (OUT, ">markers_$locilist/merlin.ped")
	or die ("open 'markers_$locilist/merlin.ped' failed, $!\n");
    while ($line = <IN>) {
	@fields = split (' ', $line);
	splice (@fields, 8, 1); # get rid of proband column
	splice (@fields, 4, 3); # and firstkid, nextpatsib and nextmatsib columns
	print (OUT join (' ', @fields), "\n");
    }
    close (IN);
    close (OUT);
}

exit (0);
