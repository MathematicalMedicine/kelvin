#!/usr/bin/perl -w
use strict;
use IO::File;

#
# seq_update_avghet - Sequentially update two-point or multipoint, Kelvin-
#                     format Bayes ratio files. 
#
# John Burian - john.burian@nationwidechildrens.org
#
# Copyright 2008, The Research Institute at Nationwide Children's Hospital
# All rights reserved. Permission is granted to use this software for
# non-profit educational purposes only.
#

my $current = "0.36.1";

my ($arg, $flag, $opt);
my %options = (relax => 0, mode => 'twopoint');
my @files;
my @filenos;
my $href;
my $key;
my $br;
my $ppl;
my $va;

my $usage = "$0: [-m {twopoint|multipoint}] [--relax] avgfile1 avgfile2 [avgfile3 ...]\n";

while ($arg = shift (@ARGV)) {
    if ($arg =~ /^--relax$/) {
	$options{relax} = 1;
    } elsif ($arg =~ /^--twopoint$/) {
        $options{mode} = 'twopoint';
    } elsif ($arg =~ /^--multipoint$/) {
        $options{mode} = 'multipoint';
    } elsif (($flag, $opt) = ($arg =~ /^-([m])(\S+)?$/)) {
        ($opt || ($opt = shift (@ARGV))) or die ($usage);
        if ($flag eq 'm') {
	    ($opt =~ /^(twopoint|multipoint)$/) or die ($usage);
            $options{mode} = $opt;
        }
    } elsif (-f $arg) {
	push (@files, { name => $arg });
    } elsif ($arg =~ /(\?|help)/) {
        system ("/usr/bin/pod2text $0");
        exit (0);
    } else {
        die ($usage);
    }
}

(scalar (@files < 2))
    and die ($usage);

foreach $href (@files) {
    ($$href{fp} = IO::File->new)->open ($$href{name})
	or die ("open '$$href{name}' failed, $!\n");
    $$href{version} = 'unknown';
    $$href{lineno} = 0;
    $$href{linetype} = '';
    $$href{expect} = 'version';
    $$href{mrknum} = '';
    $$href{mrknames} = [];
    $$href{headers} = [];
    $$href{cmpheaders} = {};
    $$href{regex} = '';
    $$href{fmt} = '';
    $$href{flds} = [];
}

@filenos = (1 .. (scalar (@files) - 1));

while (1) {
    foreach $href (@files) {
	get_next_line ($href);
	if (index ($$href{expect}, $$href{linetype}) == -1) {
	    die ("file '$$href{name}', line $$href{lineno}, expected $$href{expect}, ".
		 "found $$href{linetype}\n");
	}
    }
    
    foreach $href (@files[@filenos]) {
	($files[0]{linetype} eq $$href{linetype})
	    or die ("file '$$href{name}', line $$href{lineno}, expected $files[0]{linetype}, ".
		    "found $$href{linetype}\n");
	
	if ($$href{linetype} eq 'version') {
	    (versionnum ($$href{version}) >= versionnum ($current))
		or die ("file '$$href{name} is old format, please convert to at least $current\n");
	    $$href{expect} = 'marker or header';
	    
	} elsif ($$href{linetype} eq 'marker') {
	    unless ($options{relax} == 1) {
		($files[0]{mrknum} == $$href{mrknum}) 
		    or die ("file '$$href{name}', line $$href{lineno}, loci number mismatch\n");
		(join (',', @{$files[0]{mrknames}}) eq join (',', @{$$href{mrknames}}))
		    or die ("file '$$href{name}', line $$href{lineno}, loci name mismatch\n");
	    }
	    $$href{expect} = 'header';
	    
	} elsif ($$href{linetype} eq 'header') {
	    (join (',', sort (keys (%{$files[0]{cmpheaders}}))) eq 
	     join (',', sort (keys (%{$$href{cmpheaders}}))))
		or die ("file '$$href{name}', line $$href{lineno}, header name mismatch\n");
	    $$href{expect} = 'data';
	    #print ("headers: ", join (', ', @{$$href{headers}}), "\n");
	    #print ("cmpheaders: ", join (', ', keys (%{$$href{cmpheaders}})), "\n");
	    #print ("regex: $$href{regex}\n");
	    #print ("fmt: $$href{fmt}\n");
	    
	} elsif ($$href{linetype} eq 'data') {
	    foreach $key (keys (%{$files[0]{cmpheaders}})) {
		(($key eq 'BayesRatio') || ($key eq 'PPL')) and next;
		(($options{relax} == 1) && ($key eq 'Position')) and next;
		$va = $files[0]{cmpheaders}{$key};
		(same_value ($files[0]{flds}[$va], $$href{flds}[$va]))
		    or die ("file '$$href{name}', line $$href{lineno}, ".
			    "expected $files[0]{flds}[$va] in field '$key'\n");
	    }
	    $va = $files[0]{cmpheaders}{BayesRatio};
	    $files[0]{flds}[$va] *= $$href{flds}[$va];
	    $$href{expect} = (exists ($$href{mrknum})) ?
		'data or marker or end_of_file' : 'data or header or end_of_file';
	}
    }


    if ($files[0]{linetype} eq 'version') {
	print ("# Version V$files[0]{version}\n");
	$files[0]{expect} = 'marker or header';

    } elsif ($files[0]{linetype} eq 'marker') {
	print ("# $files[0]{mrknum} ", join (' ', @{$files[0]{mrknames}}), "\n");
	$files[0]{expect} = 'header';

    } elsif ($files[0]{linetype} eq 'header') {
	print (join (' ', @{$files[0]{headers}}), "\n");
	$files[0]{expect} = 'data';
	    
    } elsif ($files[0]{linetype} eq 'data') {
	if ($options{mode} eq 'multipoint') {
	    $va = $files[0]{cmpheaders}{BayesRatio};
	    ((($br = $files[0]{flds}[$va]) < 0.214) || 
	     (($ppl = ($br * $br) / (-5.77 + (54 * $br) + ($br * $br))) < 0))
		and $ppl = 0;
	    $va = $files[0]{cmpheaders}{PPL};
	    $files[0]{flds}[$va] = $ppl;
	}
	printf ($files[0]{fmt}, @{$files[0]{flds}});
	$files[0]{expect} = ($files[0]{mrknum}) ?
	    'data or marker or end_of_file' : 'data or header or end_of_file';
	
    } elsif ($files[0]{linetype} eq 'end_of_file') {
	exit (0);
    }
}


sub get_next_line
{
    my ($f) = @_;
    my $buff;
    my ($mrknum, $mrknames);
    my ($maj, $min, $patch, $version);
    my (@headers, $base, $numexts, @exts, $ext, $str, @regex, $keepers);
    my (@flds);
    my $va;

    my $realregex = '[\-\d\.]+(?:[eE][\+\-]\d+)?';

    $buff = $$f{fp}->getline;
    if (! defined ($buff)) {
	$$f{linetype} = 'end_of_file';
	return (1);
    }
    $$f{lineno}++;
    # These should not be necessary with versioned (non-obsolete) format BR files
    # $buff =~ s/^\s*(\S|\S.*\S)\s*$/$1/s;
    # $buff =~ s/, +/,/;

    if (($version) = ($buff =~ /\#\s+Version\s+(\S+)/)) {
	(($maj, $min, $patch) = ($version =~ /(\d+)\.(\d+)\.(\d+)/))
	    or die ("file '$$f{name}', line $$f{lineno}, malformed version number '$version'\n");
	$$f{linetype} = 'version';
	$$f{version} = sprintf ("%d.%d.%d", $maj, $min, $patch);

    } elsif (($mrknum, $mrknames) = ($buff =~ /\#\s+(\d+)((?:\s+\S+)+)/)) {
	$$f{linetype} = 'marker';
	$mrknames =~ s/^\s*(\S|\S.*\S)\s*$/$1/s;
	$$f{mrknum} = $mrknum;
	@{$$f{mrknames}} = split (/\s+/, $mrknames);

    } elsif ($buff =~ /BayesRatio/i) {
	$$f{linetype} = 'header';

	# If the regex field is non-empty, we've seen a header line
	# before, and we save ourselves a lot of work by returning now.
	# This assumes that header line (and consequently data line) 
	# formats won't change in mid-file.

	($$f{regex} ne '')
	    and return (1);
	($$f{version} ne 'unknown')
	    or die ("file '$$f{name}' is old format, please convert to at least $current\n");
	@headers = split (/\s+/, $buff);
	$va = $keepers = 0;
	@regex = ();
	while ($va < scalar (@headers)) {
	    @exts = ();

	    # Look for header fields that contain parenthesis. If a field contains
	    # an open-paren but no close-paren, tack on the next field(s), until 
	    # we find the close-paren.

	    while ($headers[$va] =~ /\([^\)]+$/) {
		($va >= scalar (@headers))
		    and die ("file '$$f{name}', line $$f{lineno}, malformed header line\n");
		$headers[$va] .= $headers[$va + 1];
		splice (@headers, $va + 1, 1);
	    }
	    $headers[$va] =~ s/\s//;

	    if (($base, $ext) = ($headers[$va] =~ /(\S+)\((\S+)*\)/)) {
		# If a header field contains parentheses, there should be a comma-
		# separated list inside that corresponds to the comma-separated
		# real numbers that will be present in the corresponding parenthesized
		# field in the data line. Count how many we should look for.
		$numexts = scalar (@exts = split (/,/, $ext));

	    } else {
		# Otherwise it's just a simple field.
		$base = $headers[$va];
		$numexts = 0;
	    }

	    if ($base eq 'Chr') {
		$$f{regex} .= '(\d+) ';
		$$f{fmt} .= '%s ';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;

	    } elsif ($base =~ /^D\d\d$/) {
		# D-Prime fields we need to keep.
		$$f{regex} .= '([\-\d\.]+) ';
		$$f{fmt} .= '%.2f ';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;
		
	    } elsif ($base eq 'Position') {
		# Position fields we need to keep.
		$$f{regex} .= '([\-\d\.]+) ';
		$$f{fmt} .= '%.4f ';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;
		
	    } elsif ($base eq 'Theta') {
		# Theta fields we need to keep.
		if ($numexts > 0) {
		    $$f{regex} .= '(\(' . '[\-\d\.]+' . ',[\-\d\.]+' x ($numexts-1) . '\)) ';
		} else {
		    $$f{regex} .= '([\-\d\.]+) ';
		}
		$$f{fmt} .= '%s ';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;
		
	    } elsif ($base eq 'PPL') {
		# PPL we match as a placeholder, and to fix the output format
		$$f{regex} .= '([\-\d\.]+) ';
		$$f{fmt} .= '%.2f ';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;
		
	    } elsif ($base eq 'BayesRatio') {
		# BayesRatio we need to keep also, but we need a more complicated regex
		$$f{regex} .= '([\-\d\.]+(?:[eE][\+\-]\d+)?) ';
		$$f{fmt} .= '%.6e ';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;
		
	    } elsif ($base =~ /(LC\d+PV|MOD|R2|Alpha|DGF|MF|MarkerList)/) {
		if ($numexts > 0) {
		    $$f{regex} .= '\(' . '[\-\d\.]+' . ',[\-\d\.]+' x ($numexts-1) . '\) ';
		    $$f{fmt} .= '(0' . ',0' x ($numexts-1) . ') ';
		} else {
		    $$f{regex} .= '[\-\d\.]+ ';
		    $$f{fmt} .= '0 ';
		}
		
	    } else {
		die ("file '$$f{name}', unrecognized column header '$base'\n");
	    }
	    push (@{$$f{headers}}, $headers[$va]);
	    $va++;
	}
	chop ($$f{fmt}, $$f{regex});
	$$f{fmt} .= "\n";

    } elsif ($$f{regex}) {
	(@flds = ($buff =~ /^$$f{regex}/))
	    or die ("file '$$f{name}', line $$f{lineno}, malformed data line\n");
	map { $flds[$_] =~ s/,\s+/,/ } (0 .. scalar (@flds)-1);
	@{$$f{flds}} = @flds;
	$$f{linetype} = 'data';
    } else {
	die ("file '$$f{name}', line $$f{lineno}, can't identify line type\n");
    }
    return (1);
}


sub same_value
{
    my ($a, $b) = @_;
    my (@a, @b);

    $a =~ s/[\(\)\s]//g;
    $b =~ s/[\(\)\s]//g;
    @a = split (/,/, $a);
    @b = split (/,/, $b);
    (scalar (@a) == scalar (@b)) or return (undef);
    while (($a = shift (@a)) && ($b = shift (@b))) {
	if ($a =~ /[\-\d\.]+(?:[Ee][\+\-]\d+)?/) {
	    ($b =~ /[\-\d\.]+(?:[Ee][\+\-]\d+)?/) or return (undef);
	    ($a == $b) or return (undef);
	} elsif ($b =~ /[\-\d\.]+(?:[Ee][\+\-]\d+)/) {
	    return (undef);
	} else {
	    ($a eq $b) or return (undef);
	}
    }
    return (1);
}


sub versionnum 
{
    my ($str) = @_;
    my ($maj, $min, $patch);

    (($maj, $min, $patch) = ($str =~ /(\d+)\.(\d+)\.(\d+)/))
	or return (0);
    return ($maj * 100000 + $min * 1000 + $patch);
}


=head1 USAGE

seq_update_br.pl [-m {twopoint|multipoint}|--twopoint|--multipoint] [--relax] brfile1 brfile2 [brfile3 ...]

Performs a sequential update of the Bayes ratios (BR)
across multilple files.

=head1 COMMAND LINE ARGUMENTS

=item -m {twopoint | multipoint}

=over 3

Sets the mode. The default is 'twopoint'.

=back

=item --twopoint

=over 3

A synonym for '-m twopoint'.

=back

=item --multipoint

=over 3

A synonym for '-m multipoint'.

=back

=item --relax

=over 3

Specifies that the loci names in the marker lines need not be the same
in all input files.

=back

=cut
