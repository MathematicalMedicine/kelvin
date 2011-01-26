#!/usr/bin/perl -w
use strict;
use KelvinIO;
use KelvinFamily 1.2;

#
# KelvinDataset: an object for managing a Kelvin-compatible marker files
# (marker map, frequencies, and locus files, and pedigree files).
#
package KelvinDataset;
our $errstr='';
our $VERSION=1.2;

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
    #   chromosome: the chromosome number from the map file
    #   mapfunction: 'kosambi' or 'haldane'
    #   writing: boolean, is the pedfile open for writing?
    #   consistent: loci data is consistent with source files
    
    #   mapfile: the name of the mapfile
    #   freqfile: the name of the allele frequency file
    #   locusfile: the name of the locus (data) file
    #   pedigreefile: the name of the pedigree file
    #   pedfh: IO::File::Kelvin handle for pedfile
    #   pedlineno: current line number in pedfile
    #   mapread: boolean, has a mapfile been read?
    #   freqread: boolean, has an allele frequency file been read?
    #   locusread: boolean, has a locus file been read?
    #   freqset: boolean, have frequencies been explicitly set?
    #   microsats: dataset contains microsatellites?
    #   snps: dataset contains snps?
    #   individualcache: cache for KelvinIndividuals when reading families
    #   pedwriteformat: 'pre' or 'post'
    
    @$self{qw/markers traits maporder/} = ({}, {}, []);
    @$self{qw/markerorder traitorder mapfields/} = ([], [], []);
    @$self{qw/chromosome mapfunction writing consistent/} = (undef, 'kosambi', 0, 1);
    @$self{qw/mapfile freqfile locusfile/} = (undef, undef, undef);
    @$self{qw/pedigreefile pedfh pedlineno/} = (undef, undef, 0);
    @$self{qw/mapread freqread locusread freqset/} = (0, 0, 0, 0);
    @$self{qw/microsats snps individualcache pedwriteformat/} = (0, 0, undef, 'post');
 
    if (! defined ($arg)) {
	$errstr = "missing argument";
	return (undef);
    } elsif (ref ($arg) eq "HASH") {
	# A hashref is allowed to specify mapfile, locusfile and freqfile
	foreach (keys (%$arg)) {
	    if (/^map/i) { $$self{mapfile} = $$arg{$_}; }
	    elsif (/^locus/i) { $$self{locusfile} = $$arg{$_}; }
	    elsif (/^freq/i) { $$self{freqfile} = $$arg{$_}; }
	    elsif (/^ped/i) { $$self{pedigreefile} = $$arg{$_}; }
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
    if (defined ($$self{pedigreefile})) {
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
    @$new{qw/chromosome mapfunction writing consistent/} = (undef, 'kosambi', 0, 1);
    @$new{qw/mapfile freqfile locusfile/} = (undef, undef, undef);
    @$new{qw/pedigreefile pedfh pedlineno/} = (undef, undef, 0);
    @$new{qw/mapread freqread locusread freqset/} = (0, 0, 0, 0);
    @$new{qw/individualcache pedwriteformat/} = (undef, 'post');

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
	if (scalar (keys (%hash)) == 0) {
	    if (scalar (@{$$self{maporder}})) {
		map { $hash{$_} = 1; } @{$$self{maporder}};
	    } else {
		map { $hash{$_} = 1; } @{$$self{markerorder}};
	    }
	}
	if (exists ($$arg{purge})) {
	    if (ref ($$arg{purge}) ne 'ARRAY') {
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
    } else {
	# Copy everything
	($$self{mapread}) and @{$$new{maporder}} = @{$$self{maporder}};
	($$self{locusread}) and @{$$new{markerorder}} = @{$$self{markerorder}};
    }
    # If markers have been dropped, mark the dataset as inconsistent with the input files
    (scalar (@{$$self{maporder}}) == scalar (@{$$new{maporder}}) &&
     scalar (@{$$self{markerorder}}) == scalar (@{$$new{markerorder}}))
	or $$new{consistent} = 0;

    # These fields always copy over, regardless of subsetting
    map {
	$$new{$_} = $$self{$_};
    } qw/chromosome mapfunction mapfile freqfile locusfile pedigreefile mapread freqread
	locusread freqset microsats snps/;

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

    } elsif ($$self{freqread} || $$self{freqset}) {
	$$new{microsats} = 0;
	$$new{snps} = 0;
	map {
	    if (exists ($$self{markers}{$_}{alleles})) {
		$href = $$self{markers}{$_}{alleles};
		@{$$new{markers}{$_}{alleles}}{keys (%$href)} = 
		    @{$$self{markers}{$_}{alleles}}{keys (%$href)};
		if (scalar (keys (%$href)) > 2) {
		    $$new{microsats} = 1;
		} else {
		    $$new{snps} = 1;
		}
	    }
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
    if ($$self{freqset}) {
	$errstr = "allele frequency file have already been set";
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

    $$self{microsats} = 0;
    $$self{snps} = 0;
    while ($line = $fh->getline (\$lineno)) {
	if ($line =~ /^M\s+(\S+)/) {
	    if ($marker) {
		if (! $self->validateAlleleFreqs ($marker, $markers{$marker}{alleles})) {
		    $errstr = "$$self{freqfile}, line ". ($lineno-1) .": $errstr";
		    return (undef);
		}
		if (scalar (keys (%{$markers{$marker}{alleles}})) > 2) {
		    $$self{microsats} = 1;
		} else {
		    $$self{snps} = 1;
		}
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
	if (scalar (keys (%{$markers{$marker}{alleles}})) > 2) {
	    $$self{microsats} = 1;
	} else {
	    $$self{snps} = 1;
	}
    } else { 
	$errstr = "$$self{freqfile} contains no markers";
	return (undef);
    }

    # If no map has been read, use the order of the markers in the frequency file
    # as a proxy for map order, so at least we can write a consistent frequency file.
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
	$errstr = "open '$$self{mapfile}' failed, $!";
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
	} elsif ($arr[$va] =~ /(sex|avg|ave|pos|kosambi)/i) {
	    push (@headers, "avgpos");
	    $regex .= (($regex) ? '\s+' : '') . '([\-\d\.]+)';
	} elsif ($arr[$va] =~ /haldane/i) {
	    push (@headers, "avgpos");
	    $$self{mapfunction} = 'haldane';
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
    $$self{chromosome} = $origchr;
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
	$$self{pedigreefile} = $pedfile;
    }
    unless (defined ($$self{pedigreefile})) {
	$errstr = "no pedigree file provided";
	return (undef);
    }
    unless ($$self{pedfh} = IO::File::Kelvin->new ($$self{pedigreefile})) {
	$errstr = "open '$$self{pedigreefile}' failed, $!";
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
    
    if ($$self{writing}) {
	$errstr = "can't read when pedigree file is open for writing";
	return (undef);
    }
    if (! $$self{consistent}) {
	$errstr = "dataset object is inconsistent with input files";
	return (undef);
    }
    (defined ($$self{pedfh}) || (defined ($$self{pedigreefile}) && $self->readPedigreefile))
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
    my $aref = [];
    my $family;
    my $pedid;
    my $pos;
    
    if (defined ($individual = $$self{individualcache})) {
	$$self{individualcache} = undef;
    } else {
	unless ($individual = $self->readIndividual) {
	    (defined ($individual)) or return (undef);
	    return (0);
	}
    }
    push (@$aref, $individual);
    $pedid = $individual->pedid;

    while ($individual = $self->readIndividual) {
	if ($individual->pedid != $pedid) {
	    $$self{individualcache} = $individual;
	    last;
	}
	push (@$aref, $individual);
    }
    (defined ($individual)) or return (undef);

    ($family = KelvinFamily->new ($aref))
	or $errstr = $KelvinFamily::errstr;
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
    $$self{consistent} = 0;
    return (1);
}

sub write
{
    my ($self, $arg) = @_;
    my ($mapfile, $locusfile, $freqfile, $pedfile);
    my $backupfile = 1;
    my $only = 1;

    if (defined ($arg)) {
	if (ref ($arg) ne "HASH") {
	    $errstr = "illegal argument";
	    return (undef);
	}
	foreach (keys (%$arg)) {
	    if (/^map/i) { $mapfile = $$arg{$_}; }
	    elsif (/^locus/i) { $locusfile = $$arg{$_}; }
	    elsif (/^freq/i) { $freqfile = $$arg{$_}; }
	    elsif (/^ped/i) { $pedfile = $$arg{$_}; }
	    elsif (/^backup$/i) { $backupfile = $$arg{$_}; }
	    elsif (/^write$/i) { 
		$only = ($$arg{$_} =~ /all/i) ? 0 : (($$arg{$_} =~ /only/i) ? 1 : undef);
	    } else {
		$errstr = "illegal argument '$_'";
		return (undef);
	    }
	}
    }
    if (! $only) {
	defined ($mapfile) or $mapfile = $$self{mapfile};
	defined ($locusfile) or $locusfile = $$self{locusfile};
        defined ($freqfile) or $freqfile = $$self{freqfile};
	defined ($pedfile) or $pedfile = $$self{pedigreefile};
    }
    if (defined ($mapfile) && $$self{mapread}) {
	$self->writeMapfile ({mapfile => $mapfile, backupfile => $backupfile})
	    or return (undef);
    }
    if (defined ($freqfile) && ($$self{freqread} || $$self{freqset})) {
	$self->writeFreqfile ({freqfile => $freqfile, backupfile => $backupfile})
	    or return (undef);
    }
    if (defined ($locusfile) && $$self{locusread}) {
	$self->writeLocusfile ({locusfile => $locusfile, backupfile => $backupfile})
	    or return (undef);
    }
    if (defined ($pedfile)) {
	$self->writePedigreefile ({pedigreefile => $pedfile, backupfile => $backupfile})
	    or return (undef);
    }
}

sub writeMapfile
{
    my ($self, $arg) = @_;
    my $mapfile = $$self{mapfile};
    my $backupfile = 1;
    my $marker;
    my $va;
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
    if ($backupfile && -f $mapfile) {
	$va = 1;
	while (-f "$mapfile.$va") {
	    $va++;
	}
	if (! rename ($mapfile, "$mapfile.$va")) {
	    $errstr = "rename '$mapfile' failed, $!";
	    return (undef);
	}
    }
    unless (open (FH, ">$mapfile")) {
	$errstr = "open '$mapfile' failed, $!";
	return (undef);
    }
    $$self{mapfile} = $mapfile;

    ($$self{mapfunction} eq 'haldane')
	and print (FH "mapFunction=haldane\n");
    print (FH join (' ', map { $headers{$_} } @{$$self{mapfields}}), "\n");
    foreach $marker (@{$$self{maporder}}) {
	print (FH join (' ', map { $$self{markers}{$marker}{$_} } @{$$self{mapfields}}), "\n");
    }
    close (FH);
    return (1);
}

sub writeFreqfile
{
    my ($self, $arg) = @_;
    my $freqfile = $$self{freqfile};
    my $backupfile = 1;
    my $marker;
    my $order;
    my $va;
    
    ($$self{freqread} || $$self{freqset}) or return (undef);
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
    if ($backupfile && -f $freqfile) {
	$va = 1;
	while (-f "$freqfile.$va") {
	    $va++;
	}
	if (! rename ($freqfile, "$freqfile.$va")) {
	    $errstr = "rename '$freqfile' failed, $!";
	    return (undef);
	}
    }
    unless (open (FH, ">$freqfile")) {
	$errstr = "open '$freqfile' failed, $!";
	return (undef);
    }
    $$self{freqfile} = $freqfile;

    $order = (scalar (@{$$self{maporder}}) > 0) ? $$self{maporder} : $$self{markerorder};
    foreach $marker (@$order) {
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
    return (1);
}

sub writeLocusfile
{
    my ($self, $arg) = @_;
    my $locusfile = $$self{locusfile};
    my $backupfile = 1;
    my $marker;
    my $trait;
    my $va;
    
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
    if ($backupfile && -f $locusfile) {
	$va = 1;
	while (-f "$locusfile.$va") {
	    $va++;
	}
	if (! rename ($locusfile, "$locusfile.$va")) {
	    $errstr = "rename '$locusfile' failed, $!";
	    return (undef);
	}
    }
    unless (open (FH, ">$locusfile")) {
	$errstr = "open '$locusfile' failed, $!";
	return (undef);
    }
    $$self{locusfile} = $locusfile;
    map { print (FH "$$self{traits}{$_}{flag} $_\n") } @{$$self{traitorder}};
    map { print (FH "M $_\n") } @{$$self{markerorder}};
    close (FH);
    return (1);
}

# This doesn't actually write the pedigree file, rather it opens the
# ped file for writing. Pedigree data is written through the 
# KelvinIndividual class.
sub writePedigreefile
{
    my ($self, $arg) = @_;
    my $pedfile = $$self{pedigreefile};
    my $backupfile = 1;
    my $premakeped = 0;
    my $marker;
    my $trait;
    my $va;

    if (defined ($arg)) {
	if (ref ($arg) eq "HASH") {
	    # A hashref is allowed to specify pedfile and backupfile
	    foreach (keys (%$arg)) {
		if (/^ped/i) { $pedfile = $$arg{$_}; }
		elsif (/^backup/i) { $backupfile = $$arg{$_}; }
		elsif (/^premakeped/i) { $premakeped = $$arg{$_}; }
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
    if ($backupfile && -f $pedfile) {
	$va = 1;
	while (-f "$pedfile.$va") {
	    $va++;
	}
	if (! rename ($pedfile, "$pedfile.$va")) {
	    $errstr = "rename '$pedfile' failed, $!";
	    return (undef);
	}
    }
    if ($backupfile && -f $pedfile && ! rename ($pedfile, "$pedfile.old")) {
	$errstr = "rename '$pedfile' failed, $!";
	return (undef);
    }
    unless ($$self{pedfh} = IO::File::Kelvin->new (">$pedfile")) {
	$errstr = "open '$pedfile' failed, $!";
	return (undef);
    }
    $$self{pedigreefile} = $pedfile;
    $$self{writing} = 1;
    $$self{pedlineno} = 0;
    ($premakeped) and $$self{pedwriteformat} = 'pre';
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
    $$self{writing} = 0;
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
    if (abs (1 - $total) <= $ROUNDING_ERROR) {
	if ($count == 1) {
	    # Assume here that a single allele, labeled either '1' or '2', with
	    # a frequency of 1, is a SNP and fill the missing complementary allele.
	    if (($label = (keys (%$href))[0]) eq '1') {
		$$href{2} = 0;
	    } elsif ($label eq '2') {
		$$href{1} = 0;
	    }
	}
    } elsif ($total < 1 - $ROUNDING_ERROR) {
	if ($count == 2) {
	    $errstr = "biallelic marker $marker frequencies don't sum to 1";
	    return (undef);
	}
	# This stinks: we really shouldn't be emitting output directly here.
	# S'okay, tho, because frequency checking will eventually move to the C code.
	warn ("WARNING, marker $marker frequencies sum to $total\n");
    } else {
	# Frquencies sum to greater than 1, and that's just wrong.
	$errstr = "marker $marker frequencies sum to $total";
	return (undef);
    }
    return (1);
}

sub setAlleleFreqs
{
    my ($self, $marker, $freqs) = @_;

    unless (exists ($$self{markers}{$marker})) {
	$errstr = "can't set allele frquencies for marker '$marker', not present in dataset";
	return (undef);
    }
    $$self{freqset} = 1;
    $$self{markers}{$marker}{alleles} = {};
    @{$$self{markers}{$marker}{alleles}}{keys (%$freqs)} = @$freqs{keys (%$freqs)};
    return (1);
}

sub mapFields
{
    my ($self) = @_;
    my $aref = [];

    @$aref = @{$$self{mapfields}};
    return ($aref);
}

sub mapOrder 
{
    my ($self) = @_;
    my $aref = [];

    @$aref = @{$$self{maporder}};
    return ($aref);
}

sub markerOrder 
{
    my ($self) = @_;
    my $aref = [];

    @$aref = @{$$self{markerorder}};
    return ($aref);
}

sub traitOrder 
{
    my ($self) = @_;
    my $aref = [];

    @$aref = @{$$self{traitorder}};
    return ($aref);
}

sub getTrait
{
    my ($self, $trait) = @_;
    my $href;

    unless (exists ($$self{traits}{$trait})) {
	$errstr = "no trait '$trait' in dataset";
	return (undef);
    }
    map { $$href{$_} = $$self{traits}{$trait}{$_} } keys (%{$$self{traits}{$trait}});
    return ($href);
}

sub pedigreefile
{
    my ($self) = @_;

    return ($$self{pedigreefile});
}

sub chromosome
{
    my ($self) = @_;

    return ($$self{chromosome});
}

sub microsats
{
    my ($self) = @_;

    return ($$self{microsats});
}

sub set_microsats
{
    my ($self, $arg) = @_;
    my $flag = ($arg) ? 1 : 0;

    return ($$self{microsats} = $flag);
}

sub snps
{
    my ($self) = @_;

    return ($$self{snps});
}

1;
