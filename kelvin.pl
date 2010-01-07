#!perl -w
use strict;
use File::Basename;
use Data::Dumper;

# KELVIN driver script 
# Copyright 2009, Nationwide Children's Hospital Research Institute
# All rights reserved. Permission is hereby granted to use this software
# for non-profit educational purposes only.

my $KELVIN_ROOT='NO_KELVIN_ROOT';
my $KELVIN_BINARY='NO_KELVIN_BINARY';
my $usage = "usage: $0 <configfile> [--directive ... ]\n";
my $config;
my $configFile;
my ($directive, $args);
my $arg;
my $idx = 0;
my $debug = 1;

if ($KELVIN_ROOT =~ /no_kelvin_root/i) {
    $KELVIN_ROOT = dirname ($0);
    print ("WARN: no KELVIN_ROOT defined, using '$KELVIN_ROOT'\n");
}
if ($KELVIN_BINARY =~ /no_kelvin_binary/i) {
    if (-x ($KELVIN_BINARY = "$KELVIN_ROOT/kelvin.". platform_name ())) {
	print ("WARN: no KELVIN_BINARY defined, using '$KELVIN_BINARY'\n");
    } elsif (-x ($KELVIN_BINARY = "$KELVIN_ROOT/kelvin.bin")) {
	print ("WARN: no KELVIN_BINARY defined, using '$KELVIN_BINARY'\n");
    } else {
	die ("FATAL: no kelvin binary found, quitting\n");
    }
}

($configFile = shift (@ARGV))
    or die ($usage);
($config = KelvinConfig->new ($configFile))
    or die ("FATAL: new KelvinConfig failed: $KelvinConfig::errstr\n");
if (@ARGV) {
    while (defined ($directive = $ARGV[$idx++])) {
	($directive =~ s/^--//) or die ($usage);
	$args = undef;
	while (defined ($ARGV[$idx]) && $ARGV[$idx] !~ /^--/) {
	    $arg = $ARGV[$idx++];
	    $args = (defined ($args)) ? $args .= " $arg" : $arg;
	}
	($config->addDirective ($directive, $args))
	    or die ("FATAL: $KelvinConfig::errstr on command line\n");
    }
    $config->validate
	or die ("FATAL: $KelvinConfig::errstr\n");
}
($debug) and print Dumper ($config);

if ($config->isConfigured ("Epistasis")) {
    run_epistasis ($config);
    
} else {
    just_run_kelvin ($config, @ARGV);
}

exit;

sub run_epistasis
{
    my ($config) = @_;
    my $epidataset;
    my $epimarker;
    my $origdataset;
    my $individual;
    my %individuals;
    my ($directive, $args);
    my $href;
    my $aref;

    # TODO: make sure we can find calc_updated_ppl

    # Configuration validation means I don't have to check return codes here
    ($directive, $args) = $config->isConfigured ("EpistasisPedigreeFile");
    $href = {pedigreefile => $$args[0]};
    ($directive, $args) = $config->isConfigured ("EpistasisLocusFile");
    $$href{locusfile} = $$args[0];
    ($epidataset = KelvinDataset->new ($href))
	or die ("FATAL: new KelvinDataset failed, $KelvinDataset::errstr\n");
    ($directive, $args) = $config->isConfigured ("Epistasis");
    $epimarker = $$args[0];
    $epidataset = $epidataset->copy ({preserve => [$epimarker]});

    while ($individual = $epidataset->readIndividual) {
	$href = $individual->structure;
	$aref = $individual->getGenotype ($epimarker);
	$individuals{"P$$href{pedid}I$$href{indid}"} = $aref;
    }
    
    ($directive, $args) = $config->isConfigured ("PedigreeFile");
    $href = {pedigreefile => $$args[0]};
    ($directive, $args) = $config->isConfigured ("LocusFile");
    $$href{locusfile} = $$args[0];
    ($origdataset = KelvinDataset->new ($href))
    	or die ("FATAL: new KelvinDataset failed, $KelvinDataset::errstr\n");


#   read locus file, add LC column, write new locus file
#   read input pedigree by individual, setting the LC in each, write first new pedigree file
#   set locus/ped file names in config, write new config
#   run kelvin with new config
}

sub just_run_kelvin
{
    my ($config, @argv) = @_;
    my $dataset;
    my $directive;
    my $args;
    
    (($directive, $args) = $config->isConfigured ("FrequencyFile"))
	or die ("no FrequencyFile in configuration file ". $config->filename);
    $dataset = KelvinDataset->new ({freqfile => $$args[0]})
	or die ($KelvinDataset::errstr. "\n");
    $dataset->close;
    print ("exec ", join (' ', map { "'$_'" } ($KELVIN_BINARY, $config->filename, @argv)), "\n");
#    exec ($KELVIN_BINARY, $config->filename, @argv)
#	or die ("exec $KELVIN_BINARY failed, $!\n");
}

sub platform_name
{
    my $platform;

    $platform = `uname -m`;
    $platform .= '-' . `uname -s`;
    $platform =~ tr/ \n/-/d;
    return ($platform);
}

my $dataset;
my $copydataset;
my $pedigree;
my $individual;
my $family;

my $href;

(($directive, $args) = $config->isConfigured ("MapFile"))
    or die ("no mapfile in configuration\n");
$href = { mapfile => $$args[0] };
(($directive, $args) = $config->isConfigured ("FrequencyFile"))
    and $$href{freqfile} = $$args[0];
(($directive, $args) = $config->isConfigured ("LocusFile"))
    and $$href{locusfile} = $$args[0];
(($directive, $args) = $config->isConfigured ("PedigreeFile"))
    and $$href{pedigreefile} = $$args[0];

$dataset = KelvinDataset->new ($href)
    or die ("$KelvinDataset::errstr\n");
($debug) and print Dumper ($dataset);

$copydataset = $dataset->copy;

$copydataset->addTrait ("Liability", "C", 'end')
    or die ("$KelvinDataset::errstr\n");
($debug) and print Dumper ($copydataset);
$copydataset->write ({ped=>'ped.new', locus=>'locus.new', map=>'map.new', freq=>'freq.new'});

($individual = $dataset->readIndividual)
    or die ("$KelvinDataset::errstr  --  $KelvinIndividual::errstr\n");
$individual->map ($copydataset);
print Dumper ($individual);
$individual->write;

#
# KelvinIndividual: an object for managing a individuals from pedigree files
#
BEGIN {
    package KelvinIndividual;
    our $errstr='';

    sub new
    {
	my ($class, $dataset, $line) = @_;
	my $trait;
	my $marker;
	my @arr;
	my $traitcol = -1;
	my $markercol = -1;
	my $colcount;
	my $ind = { pedid => undef, indid => undef, dadid => undef, momid => undef,
		    firstchildid => undef, patsibid => undef, matsibid => undef,
		    origpedid => undef, origindid => undef, sex => undef, proband => undef,
		    traits => [], markers => [], dataset => $dataset };

	# Cut off the original pedigree and person IDs, split on whitespace
	($line =~ s/\s*Ped:\s*(\w+)\s+Per:\s*(\w+)\s*$//) and
	    @$ind{qw/origpedid origindid/} = ($1, $2);
	(@$ind{qw/pedid indid dadid momid firstchildid patsibid matsibid sex proband/}, @arr) = 
	    split (' ', $line);

	# EVIL: assuming knowledge of the KelvinDataset internal structure, but so much easier
	(defined ($trait = $$dataset{traitorder}[-1]))
	    and $traitcol = $$dataset{traits}{$trait}{col};
	(defined ($marker = $$dataset{markerorder}[-1]))
	    and $markercol = $$dataset{markers}{$marker}{col};
	$colcount = ($markercol > $traitcol) ? $markercol + 2 : $traitcol + 1;

	# Basic sanity check on number of fields
	if (scalar (@arr) < $colcount) {
	    $errstr = "$$dataset{pedfile}, line $$dataset{pedlineno}: too few columns in pedigree file";
	    return (undef);
	} elsif (scalar (@arr) > $colcount && ! $$dataset{subset}) {
	    $errstr = "$$dataset{pedfile}, line $$dataset{pedlineno}: too many columns in pedigree file";
	    return (undef);
	}

	foreach $trait (@{$$dataset{traitorder}}) {
	    $traitcol = $$dataset{traits}{$trait}{col};
	    push (@{$$ind{traits}}, $arr[$traitcol]);
	}	      
	foreach $marker (@{$$dataset{markerorder}}) {
	    $markercol = $$dataset{markers}{$marker}{col};
	    push (@{$$ind{markers}}, [ $arr[$markercol], $arr[$markercol+1] ]);
	}	      
	return (bless ($ind, $class));
    }

    sub map
    {
	my ($self, $newset) = @_;
	my $oldset = $$self{dataset};
	my @traits;
	my @markers;
	my $trait;
	my $marker;

	foreach $trait (@{$$newset{traitorder}}) {
	    if (exists ($$oldset{traits}{$trait})) {
		push (@traits,  $$self{traits}[$$oldset{traits}{$trait}{idx}]);
	    } else {
		push (@traits, 'x');
	    }
	}
	foreach $marker (@{$$newset{markerorder}}) {
	    if (exists ($$oldset{markers}{$marker})) {
		push (@markers, [ @{$$self{markers}[$$oldset{markers}{$marker}{idx}]} ]);
	    } else {
		push (@markers, [ 'x', 'x' ]);
	    }
	}
	@{$$self{traits}} = @traits;
	@{$$self{markers}} = @markers;
	$$self{dataset} = $newset;
	return (1);
    }

    sub write
    {
	my ($self) = @_;
	my $dataset = $$self{dataset};
	my $origtext = '';

	(! (defined ($$dataset{pedfh}) || $$dataset{writing}))
	    and $dataset->writePedigreefile;
	(defined ($$self{origpedid}) && defined ($$self{origindid}))
	    and $origtext = "  Ped: $$self{origpedid}  Per: $$self{origindid}";

	$$dataset{pedfh}->print (join (' ', @$self{qw/pedid indid dadid momid firstchildid patsibid matsibid sex proband/}, @{$$self{traits}}, map { "$$_[0] $$_[1]" } @{$$self{markers}}), $origtext, "\n");
	return (1);
    }

    sub setTrait
    {
	my ($self, $trait, $value) = @_;
	my $dataset = $$self{dataset};
	my $idx;

	unless (exists ($$dataset{traits}{$trait})) {
	    $errstr = "no trait $trait in dataset";
	    return (undef);
	}
	$idx = $$dataset{traits}{$trait}{idx};
	(defined ($value)) or $value = 'x';
	$$self{traits}[$idx] = $value;
	return (1);
    }

    sub getGenotype
    {
	my ($self, $marker) = @_;
	my $dataset = $$self{dataset};

	(exists ($$dataset{markers}{$marker})) or return (undef);
	$idx = $$dataset{markers}{$marker}{idx};
	return ($$self{markers}[$idx]);
    }

    sub pedid
    {
	my ($self) = @_;

	return ($$self{pedid});
    }

    sub indid
    {
	my ($self) = @_;

	return ($$self{indid});
    }

    sub structure
    {
	my ($self) = @_;
	my $href = {};

	map {
	    $$href{$_} = $$self{$_};
	} qw/pedid indid dadid momid firstchildid patsibid matsibid origpedid origindid sex proband/;
	return ($href);
    }
}


#
# KelvinDataset: an object for managing a Kelvin-compatible marker files
# (marker map, frequencies, and locus files).
#
BEGIN {
    package KelvinDataset;
    #use IO::File::Kelvin;
    our $errstr='';
    my $ROUNDING_ERROR=0.0001;

    sub new
    {
	my ($class, $arg) = @_;
	my $self = bless ({}, $class);

	# Initialize the fields: 
	#   markers: a hash of all markers in the set, indexed by name
	#   traits: a hash or all traits/liability classes, indexed by name
	#   maporder: a list of markernames, ordered according to the map file
	#   markerorder: a list of markernames, ordered according to the locus file
	#   traitorder: a list of traits/liability classes, ordered according to the locus file
	#   mapfields: a list of the fields that appeared in the map
	#   mapfunction: 'kosambi' or 'haldane'
	#   writing: boolean, is the pedfile open for writing?
	#   subset: set if dataset was created as a subset of another dataset

	#   mapfile: the name of the mapfile
	#   freqfile: the name of the allele frequency file
	#   locusfile: the name of the locus (data) file
	#   pedfile: the name of the pedigree file
	#   pedfh: IO::File::Kelvin handle for pedfile
	#   pedlineno: current line number in pedfile
	#   mapread: boolean, has a mapfile been read?
	#   freqread: boolean, has an allele frequency file been read?
	#   locusread: boolean, has a locus file been read?
	#   consistant: loci data is consistant with source files
	
	@$self{qw/markers traits maporder/} = ({}, {}, []);
	@$self{qw/markerorder traitorder mapfields/} = ([], [], []);
	@$self{qw/mapfunction writing subset/} = ('kosambi', 0, 0);
	@$self{qw/mapfile freqfile locusfile/} = (undef, undef, undef);
	@$self{qw/pedfile pedfh pedlineno/} = (undef, undef, 0);
	@$self{qw/mapread freqread locusread/} = (0, 0, 0);
	@$self{qw/consistant individualcache/} = (1, undef);
	
	if (! defined ($arg)) {
	    $errstr = "missing argument";
	    return (undef);
	} elsif (ref ($arg) eq "HASH") {
	    # A hashref is allowed to specify mapfile, locusfile and freqfile
	    foreach (keys (%$arg)) {
		if (/^map/i) { $$self{mapfile} = $$arg{$_}; }
		elsif (/^locus/i) { $$self{locusfile} = $$arg{$_}; }
		elsif (/^freq/i) { $$self{freqfile} = $$arg{$_}; }
		elsif (/^ped/i) { $$self{pedfile} = $$arg{$_}; }
		else {
		    $errstr = "illegal argument '$_'";
		    return (undef);
		}
	    }
	} else {
	    # A single arg should be mapfile only. 
	    $$self{mapfile} = $arg;
	}
	if (defined ($$self{mapfile})) {
	    $self->readMapfile or return (undef);
	}
	if (defined ($$self{freqfile})) {
	    $self->readFreqfile or return (undef);
	}
	if (defined ($$self{locusfile})) {
	    $self->readLocusfile or return (undef);
	}
	if (defined ($$self{pedfile})) {
	    $self->readPedigreefile or return (undef);
	}
	return ($self);
    }

    # Create a new KelvinDataset object by copying an existing object. Copies 
    # do not inherit in the input files, if any, of the original object. Copies 
    # may be made as subsets of the original object.
    sub copy
    {
	my ($self, $arg) = @_;
	my $href;
	my $aref;
	my %hash = ();
	my $marker;
	my $idx;
	my $new = bless ({}, ref ($self));

	@$new{qw/markers traits maporder/} = ({}, {}, []);
	@$new{qw/markerorder traitorder mapfields/} = ([], [], []);
	@$new{qw/mapfunction writing subset/} = ('kosambi', 0, 0);
	@$new{qw/mapfile freqfile locusfile/} = (undef, undef, undef);
	@$new{qw/pedfile pedfh pedlineno/} = (undef, undef, 0);
	@$new{qw/mapread freqread locusread/} = (0, 0, 0);
	$$new{individualcache} = undef;

	if (defined ($arg)) {
	    if (ref ($arg) ne 'HASH') {
		$errstr = "illegal argument";
		return (undef);
	    }
	    if (exists ($$arg{preserve})) {
		if (ref ($$arg{preserve}) ne 'ARRAY') {
		    $errstr = "subset method 'preserve' requires an array reference";
		    return (undef);
		}
		map { $hash{$_} = 1; } @{$$arg{preserve}};
	    }
	    if (exists ($$arg{purge})) {
		if (ref ($$arg{preserve}) ne 'ARRAY') {
		    $errstr = "subset method 'purge' requires an array reference";
		    return (undef);
		}
		map { delete ($hash{$_}) } @{$$arg{purge}};
	    }

	    # %hash contains the final list of markers to be copied
	    map {
		exists ($hash{$_}) and push (@{$$new{maporder}}, $_);
	    } @{$$self{maporder}};
	    map {
		exists ($hash{$_}) and push (@{$$new{markerorder}}, $_);
	    } @{$$self{markerorder}};
	    $$new{subset} = 1;
	} else {
	    # Copy everything
	    ($$self{mapread}) and @{$$new{maporder}} = @{$$self{maporder}};
	    ($$self{locusread}) and @{$$new{markerorder}} = @{$$self{markerorder}};
	}

	# These fields always copy over, regardless of subsetting
	@$new{qw/mapfunction mapfile freqfile/} = @$self{qw/mapfunction mapfile freqfile/};
	@$new{qw/locusfile pedfile mapread/} = @$self{qw/locusfile pedfile mapread/};
	@$new{qw/freqread locusread consistant/} = @$self{qw/freqread locusread consistant/};

	if ($$new{mapread}) {
	    @{$$new{mapfields}} = @{$$self{mapfields}};
	    # $$self{maporder} should be non-empty here; we'll only copy the markers that
	    # appear there, in case it's a subset of the original dataset.
	    map {
		@{$$new{markers}{$_}}{@{$$new{mapfields}}} =
		    @{$$self{markers}{$_}}{@{$$new{mapfields}}};
		if (exists ($$self{markers}{$_}{alleles})) {
		    $href = $$self{markers}{$_}{alleles};
		    @{$$new{markers}{$_}{alleles}}{keys (%$href)} = 
			@{$$self{markers}{$_}{alleles}}{keys (%$href)};
		}
	    } @{$$new{maporder}};

	} elsif ($$self{freqread}) {
	    map {
		$href = $$self{markers}{$_}{alleles};
		@{$$new{markers}{$_}{alleles}}{keys (%$href)} = 
		    @{$$self{markers}{$_}{alleles}}{keys (%$href)};
	    } keys (%{$$self{markers}});
	    
	}
	if ($$new{locusread}) {
	    @{$$new{traitorder}} = @{$$self{traitorder}};
	    $idx = 0;
	    map {
		@{$$new{traits}{$_}}{qw/flag col idx/} =
		    (@{$$self{traits}{$_}}{qw/flag col/}, $idx++);
	    } @{$$new{traitorder}};
	    $idx = 0;
	    map {
		@{$$new{markers}{$_}}{qw/col idx/} = ($$self{markers}{$_}{col}, $idx++);
	    } @{$$new{markerorder}};
	}
	return ($new);
    }

    sub readLocusfile
    {
	my ($self, $locusfile) = @_;
	my $line;
	my $lineno = 0;
	my ($flag, $name);
	my $col = 0;
	my @markerorder = ();
	my @traitorder = ();
	my %markers = ();
	my %traits = ();
	my $fh;

	if (defined ($locusfile)) {
	    if (defined ($$self{locusfile})) {
		# If we've already read a locus file, forget what we read
		@{$$self{markerorder}} = ();
	    }
	    $$self{locusfile} = $locusfile;
	}
	unless (defined ($$self{locusfile})) {
	    $errstr = "no locus file provided";
	    return (undef);
	}
	unless ($fh = IO::File::Kelvin->new ($$self{locusfile})) {
	    $errstr = "open '$$self{locusfile}' failed, $!";
	    return (undef);
	}
	while ($line = $fh->getline (\$lineno)) {
	    unless (($flag, $name) = $line =~ /^(A|C|T|M)\s+(\S+)/) {
		$errstr = "$$self{locusfile}, line $lineno: badly formatted line";
		return (undef);
	    }
	    if ($flag eq 'A' || $flag eq 'T') {
		if (exists ($traits{$name})) {
		    $errstr = "$$self{locusfile}, line $lineno: trait $name appears more than once";
		    return (undef);
		}
		$traits{$name} = { col => $col, flag => $flag, idx => scalar (@traitorder) };
		push (@traitorder, $name);
		$col++;
	    } elsif ($flag eq 'C') {
		if (exists ($traits{$name})) {
		    $errstr = "$$self{locusfile}, line $lineno: liability class $name appears more than once";
		    return (undef);
		}
		$traits{$name} = {col => $col, flag => $flag, idx => scalar (@traitorder) };
		push (@traitorder, $name);
		$col++;
	    } else {  # $flag must be 'M' here
		if ($$self{mapread} && ! exists ($$self{markers}{$name})) {
		    $errstr = "$$self{locusfile}, line $lineno: marker $name in did not appear in $$self{mapfile}";
		    return (undef);
		}
		if ($$self{freqread} && ! exists ($$self{markers}{$name}{alleles})) {
		    $errstr = "$$self{locusfile}, line $lineno: marker $name in did not appear in $$self{freqfile}";
		    return (undef);
		}
		if (exists ($markers{$name})) {
		    $errstr = "$$self{locusfile}, line $lineno: marker $name appears more than once";
		    return (undef);
		}
		$markers{$name} = {col => $col, idx => scalar (@markerorder)};
		push (@markerorder, $name);
		$col += 2;
	    }
	}
	$$self{traitorder} = [ @traitorder ];
	$$self{markerorder} = [ @markerorder ];
	map {
	    @{$$self{markers}{$_}}{qw/col idx/} = @{$markers{$_}}{qw/col idx/};
	} @markerorder;
	$$self{traits} = {};
	map { $$self{traits}{$_} = $traits{$_} } @traitorder;

	$$self{locusread} = 1;
	return (1);
    }

    sub readFreqfile
    {
	my ($self, $freqfile) = @_;
	my %markers = ();
	my @maporder = ();
	my $lineno = 0;
	my ($line, @arr);
	my $marker = '';
	my $alleleno = 0;
	my $allelename = '';
	my $freq;
	my $fh;
	
	if ($$self{freqread}) {
	    $errstr = "a frequency file has already been read";
	    return (undef);
	}
	(defined ($freqfile))
	    and $$self{freqfile} = $freqfile;

	unless (defined ($$self{freqfile})) {
	    $errstr = "no frequency file provided";
	    return (undef);
	}
	unless ($fh = IO::File::Kelvin->new ($$self{freqfile})) {
	    $errstr = "open '$$self{freqfile}' failed, $!";
	    return (undef);
	}

	while ($line = $fh->getline (\$lineno)) {
	    if ($line =~ /^M\s+(\S+)/) {
		if ($marker && ! $self->validateAlleleFreqs ($marker, $markers{$marker}{alleles})) {
		    $errstr = "$$self{freqfile}, line ". ($lineno-1) .": $errstr";
		    return (undef);
		}
		push (@maporder, $marker = $1);
		if ($$self{mapread} && ! exists ($$self{markers}{$marker})) {
		    $errstr = "$$self{freqfile}, line $lineno: marker $marker in did not appear in $$self{mapfile}";
		    return (undef);
		}
		if (exists ($markers{$marker}{alleles})) {
		    $errstr = "$$self{freqfile}, line $lineno: marker $marker appears more than once";
		    return (undef);
		}
		$markers{$marker}{alleles} = {};
		$allelename = '';
		$alleleno = 0;

	    } elsif ($line =~ /^A\s+(\S+)\s+([\d\.]+)/) {
		if ($alleleno > 0) {
		    $errstr = "$$self{freqfile}, line $lineno: illegal allele specification for marker $marker";
		    return (undef);
		}
		($allelename, $freq) = ($1, $2);
		$markers{$marker}{alleles}{$allelename} = $freq;

	    } elsif (($freq) = $line =~ /^F((?:\s+[\d\.]+)+)/) {
		if ($allelename) {
		    $errstr = "$$self{freqfile}, line $lineno: illegal allele specification for marker $marker";
		    return (undef);
		}
		@arr = split (' ', $freq);
		map { $markers{$marker}{alleles}{++$alleleno} = $_; } @arr;

	    } else {
		$errstr = "$$self{freqfile}, line $lineno: badly formatted line";
		return (undef);
	    }
	}
	$fh->close;
	if ($marker) {
	    unless ($self->validateAlleleFreqs ($marker, $markers{$marker}{alleles})) {
		$errstr = "$$self{freqfile}, line $lineno: $errstr";
		return (undef);
	    }
	} else { 
	    $errstr = "$$self{freqfile} contains no markers";
	    return (undef);
	}

	# If no map has been read, use the order of the markers in the frequency file
	# as a proxy for map order, so at least we can write a consistant frequency file.
	(! $$self{mapread})
	    and @{$$self{maporder}} = @maporder;
	if ($$self{locusread}) {
	    map {
		if (! exists ($markers{$_})) {
		    $errstr = "$$self{locusfile}: marker $_ does not appear in $$self{freqfile}";
		    return (undef);
		}
	    } @{$$self{markerorder}};
	}

	foreach $marker (keys (%markers)) {
	    $$self{markers}{$marker}{alleles} = $markers{$marker}{alleles};
	}
	$$self{freqread} = 1;
	return ($self);
    }

    sub readMapfile 
    {
	my ($self, $mapfile) = @_;
	my @headers = ();
	my %markers = ();
	my @maporder = ();
	my $regex = '';
	my $lineno = 0;
	my $mapfunction;
	my $origchr = undef;;
	my $lastpos = undef;
	my ($line, @arr);
	my $href;
	my $va;
	my $fh;

	if ($$self{mapread}) {
	    $errstr = "a mapfile has already been read";
	    return (undef);
	}
	(defined ($mapfile))
	    and $$self{mapfile} = $mapfile;
	unless (defined ($$self{mapfile})) {
	    $errstr = "no map file provided";
	    return (undef);
	}
	unless ($fh = IO::File::Kelvin->new ($$self{mapfile})) {
	    $errstr = "open '$$self{mapfile} failed, $!";
	    return (undef);
	}

	# Accept one or more comment lines, and optionally a mapFunction line
	while ($line = $fh->getline (\$lineno)) {
	    ($mapfunction) = ($line =~ /mapfunction=(\w+)/i) or last;
	    if ($mapfunction =~ /^hal/i) {
		$$self{mapfunction} = 'haldane';
	    } elsif ($mapfunction !~ /^kos/i) {
		$errstr = "$$self{mapfile}, line $lineno: unknown mapFunction '$mapfunction'";
		return (undef);
	    }
	}
	@arr = split (' ', $line);
	for ($va = 0; $va < scalar (@arr); $va++) {
	    if ($arr[$va] =~ /^chr/i) {
		push (@headers, "chr");
		$regex .= (($regex) ? '\s+' : '') . '(?:chr)?(\d+)';
	    } elsif ($arr[$va] =~ /(name|marker)/i) {
		push (@headers, "name");
		$regex .= (($regex) ? '\s+' : '') . '(\S+)';
	    } elsif ($arr[$va] =~ /female/i) {
		push (@headers, "femalepos");
		$regex .= (($regex) ? '\s+' : '') . '([\-\d\.]+)';
	    } elsif ($arr[$va] =~ /male/i) {
		push (@headers, "malepos");
		$regex .= (($regex) ? '\s+' : '') . '([\-\d\.]+)';
	    } elsif ($arr[$va] =~ /(basepair|bp|phys)/i) {
		push (@headers, "phys");
		$regex .= (($regex) ? '\s+' : '') . '(\d+)';
	    } elsif ($arr[$va] =~ /(sex|avg|ave|pos)/i) {
		push (@headers, "avgpos");
		$regex .= (($regex) ? '\s+' : '') . '([\-\d\.]+)';
	    } else {
		$errstr = "$$self{mapfile}, line $lineno: unknown column header '$arr[$va]'";
		return (undef);
	    }
	}
	$line = join (' ', @headers);
	unless ($line =~ /chr/ && $line =~ /name/ && $line =~ /avgpos/) {
	    $errstr = "$$self{mapfile}, line $lineno: missing header for chromosome, marker name and/or position column";
	    return (undef);
	}

	while ($line = $fh->getline (\$lineno)) {
	    $href = {};
	    unless (@$href{@headers} = ($line =~ /^$regex/)) {
		$errstr = "$$self{mapfile}, line $lineno: badly formatted line";
		return (undef);
	    }
	    if (defined ($origchr) && $$href{chr} != $origchr) {
		$errstr = "$$self{mapfile}, line $lineno: chromosome number changes";
		return (undef);
	    }
	    if (defined ($lastpos) && $$href{avgpos} < $lastpos) {
		$errstr = "$$self{mapfile}, line $lineno: marker $$href{name} out of centiMorgan order";
		return (undef);
	    }
	    ($origchr, $lastpos) = @$href{qw/chr avgpos/};
	    $markers{$$href{name}} = $href;
	    push (@maporder, $$href{name});
	}
        $fh->close;
	if (scalar (@maporder) == 0) {
	    $errstr = "$$self{mapfile} contains no markers";
	    return (undef);
	}

	if ($$self{locusread}) {
	    map {
		if (! exists ($markers{$_})) {
		    $errstr = "$$self{locusfile}: marker $_ does not appear in $$self{mapfile}";
		    return (undef);
		}
	    } @{$$self{markerorder}};
	}
	# This works because readFreqfile will load maporder if mapread is clear
	if ($$self{freqread}) {
	    map {
		if (! exists ($markers{$_})) {
		    $errstr = "$$self{freqfile}: marker $_ does not appear in $$self{mapfile}";
		    return (undef);
		}
	    } @{$$self{maporder}};
	}
	map { 
	    @{$$self{markers}{$_}}{keys %{$markers{$_}}} = @{$markers{$_}}{keys %{$markers{$_}}};
	} keys (%markers);
	@{$$self{maporder}} = @maporder;
	@{$$self{mapfields}} = @headers;
	$$self{mapread} = 1;
	return ($self);
    }

    sub readPedigreefile
    {
	my ($self, $pedfile) = @_;	
	
	unless ($$self{locusread}) {
	    $errstr = "no locus file provided to describe pedigree file";
	    return (undef);
	}
	if (defined ($pedfile)) {
	    if (defined ($$self{pedfh})) {
		$$self{pedfh}->close;
		$$self{pedfh} = undef;
	    }
	    $$self{pedfile} = $pedfile;
	}
	unless (defined ($$self{pedfile})) {
	    $errstr = "no pedigree file provided";
	    return (undef);
	}
	unless ($$self{pedfh} = IO::File::Kelvin->new ($$self{pedfile})) {
	    $errstr = "open '$$self{pedfile}' failed, $!";
	    return (undef);
	}
	$$self{pedlineno} = 0;
	return ($self);
    }

    sub readIndividual
    {
	my ($self) = @_;
	my $line;
	my $individual;
	
	(defined ($$self{pedfh}) || (defined ($$self{pedfile}) && $self->readPedigreefile))
	    or return (undef);
	if (defined ($individual = $$self{individualcache})) {
	    $$self{individualcache} = undef;
	    return ($individual);
	}
	($line = $$self{pedfh}->getline (\$$self{pedlineno})) or return (0);
	($individual = KelvinIndividual->new ($self, $line))
	    or $errstr = $KelvinIndividual::errstr;
	return ($individual);
    }

    sub readFamily
    {
	my ($self) = @_;
	my $line;
	my $individual;
	my $family = [];
	my $pedid;
	my $pos;
	
	(defined ($$self{pedfh}) || (defined ($$self{pedfile}) && $self->readPedigreefile))
	    or return (undef);
	if (defined ($individual = $$self{individualcache})) {
	    $$self{individualcache} = undef;
	} else {
	    ($line = $$self{pedfh}->getline (\$$self{pedlineno})) or return (0);
	    unless (defined ($individual = KelvinIndividual->new ($self, $line))) {
		$errstr = $KelvinIndividual::errstr;
		return (undef);
	    }
	}
	push (@$family, $individual);
	$pedid = $individual->pedid;
	
	while ($line = $$self{pedfh}->getline (\$$self{pedlineno})) {
	    unless (defined ($individual = KelvinIndividual->new ($self, $line))) {
		$errstr = $KelvinIndividual::errstr;
		return (undef);
	    }
	    if ($individual->pedid != $pedid) {
		$$self{individualcache} = $individual;
		last;
	    }
	    push (@$family, $individual);
	}
	return ($family);
    }

    sub addTrait
    {
	my ($self, $trait, $flag, $where) = @_;
	my $va;
	my $idx;

	if (exists ($$self{traits}{$trait})) {
	    $errstr = "trait '$trait' already exists in dataset";
	    return (undef);
	}

	if ($where eq 'start') {
	    $idx = 0;
	    unshift (@{$$self{traitorder}}, $trait);
	} elsif ($where eq 'end') {
	    $idx = scalar (@{$$self{traitorder}});
	    push (@{$$self{traitorder}}, $trait);
	} else {
	    if (! exists ($$self{traits}{$where})) {
		$errstr = "can't insert trait after $where, not in dataset";
		return (undef);
	    }
	    $idx = 0;
	    # TODO: Since I already know the index of $where, why do I search for it?
	    foreach $trait (@{$$self{traitorder}}) {
		($$self{traitorder}[$idx++] eq $where) and last;
	    }
	    splice (@{$$self{traitorder}}, $idx, 0, $trait);		
	}
	$$self{traits}{$trait} = {flag => $flag};
	while ($idx < scalar (@{$$self{traitorder}})) {
	    $$self{traits}{$$self{traitorder}[$idx]}{idx} = $idx;
	    $idx++;
	}
	$$self{consistant} = 0;
	return (1);
    }

    sub write
    {
	my ($self, $arg) = @_;
	my $backupfile = 1;

	if (defined ($arg)) {
	    if (ref ($arg) ne "HASH") {
		$errstr = "illegal argument";
		return (undef);
	    }
	    foreach (keys (%$arg)) {
		if (/^map/i) { $$self{mapfile} = $$arg{$_}; }
		elsif (/^locus/i) { $$self{locusfile} = $$arg{$_}; }
		elsif (/^freq/i) { $$self{freqfile} = $$arg{$_}; }
		elsif (/^ped/i) { $$self{pedfile} = $$arg{$_}; }
		elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		else {
		    $errstr = "illegal argument '$_'";
		    return (undef);
		}
	    }
	}
	if (defined ($$self{mapfile}) && $$self{mapread}) {
	    $self->writeMapfile ({backupfile => $backupfile}) or return (undef);
	}
	if (defined ($$self{freqfile}) && $$self{freqread}) {
	    $self->writeFreqfile ({backupfile => $backupfile}) or return (undef);
	}
	if (defined ($$self{locusfile}) && $$self{locusread}) {
	    $self->writeLocusfile ({backupfile => $backupfile}) or return (undef);
	}
	if (defined ($$self{pedfile})) {
	    $self->writePedigreefile ({backupfile => $backupfile}) or return (undef);
	}
    }

    sub writeMapfile
    {
	my ($self, $arg) = @_;
	my $mapfile = $$self{mapfile};
	my $backupfile = 1;
	my $marker;
	my %headers = (chr => 'Chromosome', name => 'Marker',
		       femalepos => 'FemalePosition', malepos => 'MalePosition',
		       avgpos => 'Position', phys => 'Basepair');
	
	($$self{mapread}) or return (undef);
	if (defined ($arg)) {
	    if (ref ($arg) eq "HASH") {
		# A hashref is allowed to specify mapfile and backupfile
		foreach (keys (%$arg)) {
		    if (/^map/i) { $mapfile = $$arg{$_}; }
		    elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		    else {
			$errstr = "illegal argument '$_'";
			return (undef);
		    }
		}
	    } else {
		# A single arg should be mapfile only. 
		$mapfile = $arg;
	    }
	}
	unless (defined ($mapfile)) {
	    $errstr = "no map file specified";
	    return (undef);
	}
	if ($backupfile && -f $mapfile && ! rename ($mapfile, "$mapfile.old")) {
	    $errstr = "rename '$mapfile' failed, $!";
	    return (undef);
	}
	unless (open (FH, ">$mapfile")) {
	    $errstr = "open '$mapfile' failed, $!";
	    return (undef);
	}
	print (FH join (' ', map { $headers{$_} } @{$$self{mapfields}}), "\n");
	foreach $marker (@{$$self{maporder}}) {
	    print (FH join (' ', map { $$self{markers}{$marker}{$_} } @{$$self{mapfields}}), "\n");
	}
	close (FH);
    }

    sub writeFreqfile
    {
	my ($self, $arg) = @_;
	my $freqfile = $$self{freqfile};
	my $backupfile = 1;
	my $marker;
	
	($$self{freqread}) or return (undef);
	if (defined ($arg)) {
	    if (ref ($arg) eq "HASH") {
		# A hashref is allowed to specify freqfile and backupfile
		foreach (keys (%$arg)) {
		    if (/^freq/i) { $freqfile = $$arg{$_}; }
		    elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		    else {
			$errstr = "illegal argument '$_'";
			return (undef);
		    }
		}
	    } else {
		# A single arg should be freqfile only. 
		$freqfile = $arg;
	    }
	}
	unless (defined ($freqfile)) {
	    $errstr = "no frequency file specified";
	    return (undef);
	}
	if ($backupfile && -f $freqfile && ! rename ($freqfile, "$freqfile.old")) {
	    $errstr = "rename '$freqfile' failed, $!";
	    return (undef);
	}
	unless (open (FH, ">$freqfile")) {
	    $errstr = "open '$freqfile' failed, $!";
	    return (undef);
	}

	foreach $marker (@{$$self{maporder}}) {
	    (exists ($$self{markers}{$marker}{alleles})) or next;
	    print (FH "M $marker\n");
	    if (join (keys (%{$$self{markers}{$marker}{alleles}})) =~ /^\d+$/) {
		map {
		    print (FH "A $_ $$self{markers}{$marker}{alleles}{$_}\n");
		} sort ({$a <=> $b} keys (%{$$self{markers}{$marker}{alleles}}));
	    } else {
		map {
		    print (FH "A $_ $$self{markers}{$marker}{alleles}{$_}\n");
		} sort (keys (%{$$self{markers}{$marker}{alleles}}));
	    }
	}
	close (FH);
    }

    sub writeLocusfile
    {
	my ($self, $arg) = @_;
	my $locusfile = $$self{locusfile};
	my $backupfile = 1;
	my $marker;
	my $trait;
	
	($$self{locusread}) or return (undef);
	if (defined ($arg)) {
	    if (ref ($arg) eq "HASH") {
		# A hashref is allowed to specify locusfile and backupfile
		foreach (keys (%$arg)) {
		    if (/^locus/i) { $locusfile = $$arg{$_}; }
		    elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		    else {
			$errstr = "illegal argument '$_'";
			return (undef);
		    }
		}
	    } else {
		# A single arg should be locusfile only. 
		$locusfile = $arg;
	    }
	}
	unless (defined ($locusfile)) {
	    $errstr = "no locus file specified";
	    return (undef);
	}
	if ($backupfile && -f $locusfile && ! rename ($locusfile, "$locusfile.old")) {
	    $errstr = "rename '$locusfile' failed, $!";
	    return (undef);
	}
	unless (open (FH, ">$locusfile")) {
	    $errstr = "open '$locusfile' failed, $!";
	    return (undef);
	}
	map { print (FH "$$self{traits}{$_}{flag} $_\n") } @{$$self{traitorder}};
	map { print (FH "M $_\n") } @{$$self{markerorder}};
	close (FH);
    }

    sub writePedigreefile
    {
	my ($self, $arg) = @_;
	my $pedfile = $$self{pedfile};
	my $backupfile = 1;
	my $marker;
	my $trait;
	
	if (defined ($arg)) {
	    if (ref ($arg) eq "HASH") {
		# A hashref is allowed to specify pedfile and backupfile
		foreach (keys (%$arg)) {
		    if (/^ped/i) { $pedfile = $$arg{$_}; }
		    elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		    else {
			$errstr = "illegal argument '$_'";
			return (undef);
		    }
		}
	    } else {
		# A single arg should be pedfile only. 
		$pedfile = $arg;
	    }
	}
	unless (defined ($pedfile)) {
	    $errstr = "no pedigree file specified";
	    return (undef);
	}
	(defined ($$self{pedfh}))
	    and $$self{pedfh}->close;
	if ($backupfile && -f $pedfile && ! rename ($pedfile, "$pedfile.old")) {
	    $errstr = "rename '$pedfile' failed, $!";
	    return (undef);
	}
	unless ($$self{pedfh} = IO::File::Kelvin->new (">$$self{pedfile}")) {
	    $errstr = "open '$$self{pedfile}' failed, $!";
	    return (undef);
	}
	$$self{writing} = 1;
	$$self{pedlineno} = 0;
	return (1);
    }

    sub close 
    {
	my ($self) = @_;

	if (defined ($$self{pedfh})) {
	    $$self{pedfh}->close;
	    $$self{pedfh} = undef;
	}
	$$self{individualcache} = undef;
	return (1);
    }

    sub validateAlleleFreqs
    {
	my ($self, $marker, $href) = @_;
	my $count = 0;
	my $total = 0;
	my $label;

	map { $count++; $total += $$href{$_} } keys %$href;
	if ($count == 0) {
	    $errstr = "marker $marker has no alleles";
	    return (undef);
	}
	if (abs (1 - $total) < $ROUNDING_ERROR) {
	    if ($count == 1) {
		if (($label = (keys (%$href))[0]) eq '1') {
		    $$href{2} = 0;
		} elsif ($label eq '2') {
		    $$href{1} = 0;
		}
	    }
	    return (1);
	}
	if ($count == 2) {
	    $errstr = "biallelic marker frequencies don't sum to 1";
	    return (undef);
	}
	return (1);
    }

    sub mapFields
    {
	my ($self) = @_;

	return ($$self{mapfields});
    }

    sub mapOrder 
    {
	my ($self) = @_;

	return ($$self{maporder});
    }

    sub markerOrder 
    {
	my ($self) = @_;

	return ($$self{markerorder});
    }

    sub traitOrder 
    {
	my ($self) = @_;

	return ($$self{traitorder});
    }

    sub markers
    {
	my ($self) = @_;

	return ($$self{markers});
    }

    sub traits
    {
	my ($self) = @_;

	return ($$self{traits});
    }
}


#
# KelvinConfig: an object for managing a Kelvin configuration file.
#
BEGIN {
    package KelvinConfig;
    #use IO::File::Kelvin;
    our $errstr='';
    my %directives = (
		      pedigreefile => {canon => 'PedigreeFile',
				       singlearg => 'true'},
		      locusfile => {canon => 'LocusFile',
				    singlearg => 'true'},
		      frequencyfile => {canon => 'FrequencyFile',
					singlearg => 'true'},
		      mapfile => {canon => 'MapFile',
				  singlearg => 'true'},
		      bayesratiofile => {canon => 'BayesRatioFile',
					 singlearg => 'true'},
		      pplfile => {canon => 'PPLFile',
				  singlearg => 'true'},
		      countfile => {canon => 'CountFile',
				    singlearg => 'true'},
		      modfile => {canon => 'MODFile',
				  singlearg => 'true'},
		      extramods => {canon => 'ExtraMODs'},
		      forcebrfile => {canon => 'ForceBRFile'},
		      surfacespath=> {canon => 'SurfacesPath',
				      singlearg => 'true'},
		      surfacefile => {canon => 'SurfaceFile',
				      singlearg => 'true'},
		      nidetailfile => {canon => 'NIDetailFile',
				       singlearg => 'true'},
		      epistasispedigreefile => {canon => 'EpistasisPedigreeFile',,
						singlearg => 'true',
						local => 'true',
						regex => '(\S+)'},
		      epistasislocusfile => {canon => 'EpistasisLocusFile',,
					     singlearg => 'true',
					     local => 'true',
					     regex => '(\S+)'},

		      multipoint => {canon => 'MultiPoint'},
		      markertomarker => {canon => 'MarkerToMarker'},
		      ld => {canon => 'LD'},
		      epistasis => {canon => 'Epistasis',
				    local => 'true',
				    regex => '([\w\-]+(?:,\s*[\w\-]+)*)'},


		      qt => {canon => 'QT'},
		      qtt => {canon => 'QTT'},
		      threshold => {canon => 'Threshold'},
		      liabilityclasses => {canon => 'LiabilityClasses'},
		      diseasegenefrequency => {canon => 'DiseaseGeneFrequency'},
		      dprime=> {canon => 'DPrime'},
		      theta => {canon => 'Theta'},
		      alpha => {canon => 'Alpha'},
		      penetrance => {canon => 'Penetrance'},
		      mean => {canon => 'Mean'},
		      standarddev => {canon => 'StandardDev'},
		      degreesoffreedom => {canon => 'DegreesOfFreedom'},
		      constraint => {canon => 'Constraint'},
		      truncate=> {canon => 'Truncate'},
		      markerallelefrequency => {canon => 'MarkerAlleleFrequency'},

		      phenocodes => {canon => 'PhenoCodes'},
		      sexspecific => {canon => 'SexSpecific'},
		      sexlinked => {canon => 'SexLinked'},
		      imprinting => {canon => 'Imprinting'},
		      traitpositions => {canon => 'TraitPositions'},
		      diseasealleles=> {canon => 'DiseaseAlleles'},
		      polynomialscale=> {canon => 'PolynomialScale'},
		      nonpolynomial=> {canon => 'NonPolynomial'},
		      fixedmodels => {canon => 'FixedModels'},
		      dryrun => {canon => 'DryRun'},
		      maxiterations => {canon => 'MaxIterations'},
		      log=> {canon => 'Log',
			     regex => '(\w+)\s+(\w+)'}
		      );
    
    sub new
    {
	my ($class, $configfile) = @_;
	my $self;
	my $line;
	my $lineno = 0;
	my ($directive, $arg);
	my $fh;

	unless ($fh = IO::File::Kelvin->new ($configfile)) {
	    $errstr = $!;
	    return (undef);
	}
	$self = bless ({ filename => $configfile, directives => {} }, $class);
	while ($line = $fh->getline (\$lineno)) {
	    while ($line =~ s/(\w+)(?:\s+([^;\n\r]+))?;?//) {
		($directive, $arg) = ($1, $2);
		unless ($directive = $self->addDirective ($directive, $arg)) {
		    $errstr .= " on line $lineno";
		    return (undef);
		}
	    }
	}
	$fh->close;
	($self->validate) or return (undef);
	return ($self);
    }

    sub addDirective
    {
	my ($self, $directive, $arg) = @_;
	
	($directive = $self->legalDirective ($directive))
	    or return (undef);
	if (exists ($directives{lc($directive)}{regex}) &&
	    (! defined ($arg) || $arg !~ /^$directives{lc($directive)}{regex}$/)){
	    $errstr = "Illegal argument to directive $directive";
	    return (undef);
	}
	exists ($$self{directives}{$directive}) or 
	    $$self{directives}{$directive} = [];
	if (defined ($arg)) {
	    if (exists ($directives{lc($directive)}{singlearg}) && 
		$directives{lc($directive)}{singlearg} eq 'true') {
		$$self{directives}{$directive}[0] = $arg;
	    } else {
		push (@{$$self{directives}{$directive}}, $arg);
	    }
	}
	return ($directive);
    }

    sub removeDirective
    {
	my ($self, $directive) = @_;

	($directive = $self->legalDirective ($directive))
	    or return (undef);
	unless (exists ($$self{directives}{$directive})) {
	    $errstr = "$directive is not in configuration";
	    return (undef);
	}
	delete ($$self{directives}{$directive});
	return (1);
    }

    sub isConfigured
    {
	my ($self, $directive) = @_;

	($directive = $self->legalDirective ($directive))
	    or return (undef);
	exists ($$self{directives}{$directive})
	    or return (undef);
	return ($directive, $$self{directives}{$directive});
    }

    sub write
    {
	my ($self, $arg) = @_;
	my $configfile = undef;
	my $backupfile = 1;
	my $nolocal = 0;
	my $directive;

	if (! defined ($arg)) {
	    $errstr = "missing argument";
	    return (undef);
	} elsif (ref ($arg) eq "HASH") {
	    # A hashref is allowed to specify configfile, backupfile and nolocal
	    foreach (keys (%$arg)) {
		if (/^config/i) { $configfile = $$arg{$_}; }
		elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		elsif (/^nolocal/i) { $nolocal = $$arg{$_}; }
		else {
		    $errstr = "illegal argument '$_'";
		    return (undef);
		}
	    }
	} else {
	    # A single arg should be configfile only. 
	    $configfile = $arg;
	}
	unless (defined ($configfile)) {
	    $errstr = "no configfile specified";
	    return (undef);
	}
	if ($backupfile && -f $configfile && ! rename ($configfile, "$configfile.old")) {
	    $errstr = "rename '$configfile' failed, $!";
	    return (undef);
	}
	unless (open (FH, ">$configfile")) {
	    $errstr = "open '$configfile' failed, $!";
	    return (undef);
	}
	foreach $directive (keys (%{$$self{directives}})) {
	    ($nolocal && exists ($directives{lc($directive)}{local}) &&
	     $directives{lc($directive)}{local} eq 'true')
		and next;
	    if (scalar (@{$$self{directives}{$directive}}) == 0) {
		print (FH $directive, "\n");
	    } else {
		map { print (FH $directive . ' ' . $_, "\n"); } @{$$self{directives}{$directive}};
	    }
	}
	close (FH);
	return (1);
    }

    sub legalDirective
    {
	my ($self, $nominal) = @_;
	my $directive;
	my @matches = ();

	unless (defined ($nominal)) {
	    $errstr = "missing argument";
	    return (undef);
	}
	$nominal = lc ($nominal);
	exists ($directives{$nominal})
	    and return ($directives{$nominal}{canon});
	foreach $directive (keys (%directives)) {
	    ($directive =~ /^$nominal/)
		and push (@matches, $directive);
	}
	if (scalar (@matches) == 0) {
	    $errstr = "illegal directive '$nominal'";
	    return (undef);
	} elsif (scalar (@matches) > 1) {
	    $errstr = "ambiguous directive '$nominal'";
	    return (undef);
	} else {
	    return ($directives{$matches[0]}{canon});
	}
    }

    sub validate
    {
	my ($self) = @_;
	my $directive;
	
	# We only care about local directives here; the Kelvin executable 
	# can handle regular directives.

	# Epistasis, EpistasisPedigreeFile and EpistatisLocusFile are mutually
	# dependent, and incompatible with LiabilityClasses
	foreach $directive qw(Epistasis EpistasisPedigreeFile EpistasisLocusFile) {
	    exists ($$self{directives}{$directive}) or next;
	    map {
		if (! exists ($$self{directives}{$_})) {
		    $errstr = "$directive requires $_";
		    return (undef);
		}
	    } qw(Epistasis EpistasisPedigreeFile EpistasisLocusFile);
	    if (exists ($$self{directives}{LiabilityClasses})) {
		$errstr = "$directive is incompatible with LiabilityClasses";
		return (undef);
	    }
	}
	
	return (1);
    }

    sub filename
    {
	my ($self) = @_;
	
	return ($$self{filename});
    }
}


#
# IO::File::Kelvin: A convenience wrapper around IO::FIle that only returns
# nonblank, noncomment lines, and counts lines read.
#
BEGIN {
    package IO::File::Kelvin;
    require IO::File;

    our @ISA = qw(IO::File);

    sub new
    {
	my ($class, @args) = @_;
	return ($class->SUPER::new(@args));
    }

    sub getline
    {
	my ($self, $ref) = @_;
        my $line;
	my $lineno=0;
	
	while ($line = $self->SUPER::getline) {
	    $lineno++;
	    $line =~ s/^\s*([^\#\n]*)\#?.*$/$1/s;
	    if ($line) {
		# This is clumsy. How to do it better?
		(ref ($ref) eq 'SCALAR' && $$ref =~ /^\d+$/) and $$ref += $lineno;
		return ($line);
	    }
        }
	return (undef);
    }
}
