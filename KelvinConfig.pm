#!/usr/bin/perl -w
use strict;
use KelvinIO;

=head1 SUMMARY

KelvinConfig provides an object for reading, querying and writing Kelvin
configuration files. Only limited configuration validation is performed, 
generally what can be accomplished by matching directive arguments 
against simple regular expressions. This object is largely written for
use in the Perl frontend for the Kelvin package. As such, it makes 
special considerations for directives that are used by the frontend.
Detection of conflicting directives is only performed for directives
that are handled directly in the frontend.

=head1 USAGE

	use KelvinConfig;
	$config = KelvinConfig->new ($filename);
	$config->addDirective ($directive, $arg);
	$config->setDirective ($directive, $arg);
	$config->removeDirective ($directive);
	$config->legalDirective ($directive);
	$config->write ($filename);
	$config->write ($href);
	$config->filename;
	$config->localDirectives;

In all cases, $directive is case-insensitive, and may be an abbreviation,
so long as it unambiguously matches exactly one legal directive. In all
cases, $arg should be a simple string, such as would appear in a Kelvin
configuration file after a directive. All methods set $errstr and return
undef on failure.

=head1 METHODS

=cut

package KelvinConfig;
our $errstr='';
our $VERSION=1.0;


# The %directive hash: The keys are lowercased versions of the 
#   canonical directive name. Each key must point to a hash ref
#   that itself must contain at least one key, 'canon', which
#   points to the canonical directive name. The hash ref may
#   optionally contain additional keys:
#
#     singlearg: If the value is 'true', only a single argument is
#       applicable to the directive. The 'addDirective' method will
#       behave like the 'setDirective' method, overwriting any
#       previous argument for the directive, rather than appending.
#     default: If present, the value is the default argument for
#       the directive. If there is no default argument, this key should
#       not exist. If the default argument is universally applicable,
#       the value can be a plain string. If the default depends on the
#       presence of, or arguments to, other direcrtives, then this
#       can be a reference to a method that will return the correct
#       default argument.
#     local: If the value is 'true', then the directive is only legal
#       in the Perl frontend. The C code will not accept such directives.
#     regex: The value is a Perl regular expression that the directive
#       argument must match. Note that the regex must be flexible enough
#       to match complex arguments, like comma-separated lists, optional
#       keywords, etc. If the regex matches, the entire argument string
#       will be stored; no consideration is given to parenthesized
#       sub-expressions.
#     parser: The value is a reference to a method that will parse the
#       argument. The method must return either a reference to a 
#       scalar, or a reference to an array.
#
# Argument regular-expression matching (the 'regex' key) or parsing
# the 'parser' key) should be used sparingly, and only in cases where
# KelvinConfig methods or the Perl frontend require it. The 'regex'
# and 'parser' keys should never both exist for a directive. In general,
# the KelvinConfig object should do as little validation as possible,
# since the compiled Kelvin binary will exhaustively validate the 
# config.

my %directives = (
		  pedigreefile => {canon => 'PedigreeFile',
				   singlearg => 'true',
				   default => ['pedfile.dat']},
		  locusfile => {canon => 'LocusFile',
				singlearg => 'true',
				default => ['datafile.dat']},
		  frequencyfile => {canon => 'FrequencyFile',
				    singlearg => 'true'},
		  mapfile => {canon => 'MapFile',
			      singlearg => 'true',
			      default => ['mapfile.dat']},
		  bayesratiofile => {canon => 'BayesRatioFile',
				     singlearg => 'true',
				     default => ['br.out']},
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
		  epistasispedigreefile => {canon => 'EpistasisPedigreeFile',
					    singlearg => 'true',
					    local => 'true',
					    regex => '\S+'},
		  epistasislocusfile => {canon => 'EpistasisLocusFile',
					 singlearg => 'true',
					 local => 'true',
					 regex => '\S+'},

		  multipoint => {canon => 'MultiPoint'},
		  markertomarker => {canon => 'MarkerToMarker',
				     singlearg => 'true',
				     default => ['adjacent'],
				     regex => '(?:all|adjacent)'},
		  ld => {canon => 'LD'},
		  epistasis => {canon => 'Epistasis',
				local => 'true',
				regex => '[\w\-]+(?:,\s*[\w\-]+)*'},

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

		  # Note we validate args to PhenoCodes in validate()
		  phenocodes => {canon => 'PhenoCodes',
				 singlearg => 'true',
				 default => \&defaultPhenoCodes,
				 regex => '[\-\d\.]+(?:\s*,\s*[\-\d\.]+){2}'},
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
			 regex => '\w+\s+\w+'},
		  skippedcount => {canon => 'SkipPedCount',
				   local => 'true'},
		  skipanalysis => {canon => 'SkipAnalysis',
				   local => 'true'},
		  generatekeywords => {canon => 'GenerateKeywords',
				       local => 'true'},
		  skipestimation => {canon => 'SkipEstimation',
					local => 'true'},
		  study => {canon => 'Study',
			    local => 'true',
			    singlearg => 'true',
			    regex => '(\d+)\s+(\w+)\s+\"(.+)\"\s+(\w+)\s+(\w+)'},
		  traitprevalence => {canon => 'TraitPrevalence',
				      local => 'true',
				      singlearg => 'true',
				      regex => '0*\.[0-9]+'}
		  );


=over 2

=item $config = KelvinConfig->new ($filename);

Creates a KelvinConfig object based on the contents of $filename, which must
be a Kelvin 2.0+ (verbose config style) configuration file. Returns the 
object reference on success. Success only guarantees that the file contained
legal directives, it does not guarantee that all arguments are correct or 
that conflicting directives or arguments were not present.

=back 

=cut

sub new
{
    my ($class, $configfile) = @_;
    my $self;
    my $line;
    my $lineno = 0;
    my ($directive, $arg);
    my $fh;

    unless (defined ($configfile)) {
	$errstr = "config file name is undefined.";
	return (undef);
    }
    unless ($fh = IO::File::Kelvin->new ($configfile)) {
	$errstr = $!;
	return (undef);
    }
    $self = bless ({ filename => $configfile, directives => {} }, $class);
    $$self{localdirectives} = 0;
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

=over 2

=item $newconfig = $config->copy;

Creates a new KelvinConfig object by copying $config.

=back 

=cut

sub copy
{
    my ($self) = @_;
    my $new = {};
    my $directive;

    map { $$new{$_} = $$self{$_} } qw/filename localdirectives/;
    $$new{directives} = {};

    foreach $directive (keys (%{$$self{directives}})) {
	$$new{directives}{$directive} = [];
	@{$$new{directives}{$directive}} = @{$$self{directives}{$directive}};
    }
    return (bless ($new, ref ($self)));
}

=over 2

=item $bool = $config->addDirective ($directive, $arg);

Adds a directive to the $config. If $directive does not already 
exists in $config, the directive will be added. If $arg is provided,
it will be added to the list of argument strings for the directive,
unless the directive is of the single-argument type, in which case
$arg will replace any previously configured arguments. Returns
the canonical version of $directive on success.

=back 

=cut

sub addDirective
{
    my ($self, $directive, $arg) = @_;
    my $ref = undef;
    
    ($directive = $self->legalDirective ($directive))
	or return (undef);
    unless (exists ($$self{directives}{$directive})) {
	$$self{directives}{$directive} = [];
	($directives{lc($directive)}{local})
	    and $$self{localdirectives}++;
    }
    if (defined ($arg)) {
	if (exists ($directives{lc($directive)}{regex}) &&
	    $arg !~ /^$directives{lc($directive)}{regex}$/i) {
	    $errstr = "Illegal argument to directive $directive";
	    return (undef);
	} elsif (exists ($directives{lc($directive)}{parser}) &&
		 ! ($ref = &{$directives{lc($directive)}{parser}} ($self, $arg))) {
	    $errstr = "Illegal argument to directive $directive";
	    return (undef);
	}
	if (exists ($directives{lc($directive)}{singlearg}) && 
	    $directives{lc($directive)}{singlearg} eq 'true') {
	    if (! defined ($ref)) {
		$$self{directives}{$directive}[0] = $arg;
	    } elsif (ref ($ref) eq 'SCALAR') {
		$$self{directives}{$directive}[0] = $$ref;
	    } else { # ref to an ARRAY
		@{$$self{directives}{$directive}} = @$ref;
	    }
	} else {
	    if (! defined ($ref)) {
		push (@{$$self{directives}{$directive}}, $arg);
	    } elsif (ref ($ref) eq 'SCALAR') {
		push (@{$$self{directives}{$directive}}, $$ref);
	    } else { # ref to an ARRAY
		push (@{$$self{directives}{$directive}}, @$ref);
	    }
	}
    } else {
	if (exists ($directives{lc($directive)}{regex}) ||
	    exists ($directives{lc($directive)}{parser})) {
	    $errstr = "Illegal argument to directive $directive";
	    return (undef);
	}
    }
    return ($directive);
}

=over 2

=item $bool = $config->setDirective ($directive, $arg);

Does the same thing as addDirective, but replaces any previous
argument(s) to $directive with the supplied $arg, regardless of
whether or not the directive allows multiple arguments. Returns
the canonical version of $directive on success.

=back 

=cut

sub setDirective
{
    my ($self, $directive, $arg) = @_;

    ($directive = $self->legalDirective ($directive))
	or return (undef);
    (exists ($$self{directives}{$directive}))
	and $$self{directives}{$directive} = [];
    return ($self->addDirective ($directive, $arg));
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
    ($directives{lc($directive)}{local})
	and $$self{localdirectives}--;
    delete ($$self{directives}{$directive});
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
    
    # For the most part, we only care about local directives here; the Kelvin
    # executable  can handle regular directives.

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

    # We have to validate PhenoCodes here because the front end requires it.
    if (exists ($$self{directives}{PhenoCodes}) && exists ($$self{directives}{QT})) {
	if ($$self{directives}{PhenoCodes}[0] !~ /^[\-\d\.]+$/) {
	    $errstr = "Illegal argument to directive PhenoCodes";
	    return (undef);
	}
    }
    return (1);
}

=over 2

=item $bool = $config->write ($filename);

=item $bool = $config->write ($hashref);

Writes the configuration to a Kelvin 2.0+ format config file. In the first 
form, the config is written to $filename, with default options (see below).
In the second form, the filename and options are specified by passing a
reference to a hash. The hash may contain the following keys:

=over 2

configfile - the value is the filename to which the config will be written.

backupfile - if the value is true, and the config file already exists, it
will be renamed with the extension '.old'.

nolocal - if the value is true, local directives (directives that are only
legal in the Perl frontend) will not be written to the new config file.

=back

By default, configfile is the file from which the config was read, backupfile
is true, and nolocal is false.

=back 

=cut

sub write
{
    my ($self, $arg) = @_;
    my $configfile = $$self{filename};
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
    $$self{filename} = $configfile;
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

=over 2

=item $string = $config->filename;

Returns the filename from which the config was read.

=back 

=cut

sub filename
{
    my ($self) = @_;
    
    return ($$self{filename});
}

=over 2

=item $bool = $config->localDirectives;

Returns true if $config contains any local directives (directives that
are only legal in the Perl frontend), false otherwise.

=back 

=cut

sub localDirectives
{
    my ($self) = @_;

    return ($$self{localdirectives});
}

=over 2

=item $arrayref = $config->isConfigured ($directive);

=item ($directive, $arrayref) = $config->isConfigured ($directive);

Determines if $directive exists in $config. In a scalar context,
returns an array reference containing the list of configured
argument strings, if any. In an array context, returns the canonical
version of $directive, and the aforementioned array reference. If
the directive is not configured, but there is a default argument,
the default will be returned in the array reference.

=back 

=cut

sub isConfigured
{
    my ($self, $directive) = @_;
    my $aref;

    ($directive = $self->legalDirective ($directive))
	or return (undef);
    if (exists ($$self{directives}{$directive})) {
	$aref = $$self{directives}{$directive};
    } elsif (exists ($directives{lc($directive)}{default})) {
	$aref = (ref ($directives{lc($directive)}{default}) eq 'CODE') ?
	    &{$directives{lc($directive)}{default}} ($self) :
	    $directives{lc($directive)}{default};
    } else {
	return (undef);
    }
    if (wantarray ()) {
	return ($directive, $aref);
    } else {
	return ($aref);
    }
}

sub defaultPhenoCodes
{
    my ($self) = @_;

    if (exists ($$self{directives}{'QT'})) {
	# QT only has a phenocode for 'unknown'
	return (['-99.99']);
    } elsif (exists ($$self{directives}{'QTT'})) {
	# QTT codes for unknown, unaffected and affected
	return (['-99.99, -88.88, 88.88']);
    } else {
	# DT codes for unknown, unaffected and affected
	return (['0, 1, 2']);
    }
}

1;