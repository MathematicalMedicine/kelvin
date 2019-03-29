#!/usr/bin/perl -w
use strict;
use IO::File;


my ($arg, $flag, $opt);
my %options = (relax => 0, mode => 'twopoint');
my @files;
my @filenos;
my $href;
my $key;
my $avglr;
my $ppl;
my $va;

my $usage = "$0: [-m mode] [--relax] avgfile1 avgfile2 [avgfile3 ...]\n";

while ($arg = shift (@ARGV)) {
    if ($arg =~ /^--relax$/) {
	$options{relax} = 1;
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
    $$href{lineno} = 0;
    $$href{linetype} = '';
    $$href{expect} = 'marker line';
    $$href{mrknum} = '';
    $$href{mrknames} = [];
    $$href{headers} = [];
    $$href{cmpheaders} = {};
    $$href{regex} = '';
    $$href{format} = '';
    $$href{flds} = [];
}

@filenos = (1 .. (scalar (@files) - 1));

while (1) {
    foreach $href (@files) {
	get_next_line ($href);
#	($$href{linetype} != $$href{expect})
#	    or die ("file '$$href{name}', line $$href{lineno}, expected $$href{expect}, ".
#		    "found $$href{linetype}\n");
    }
    
    foreach $href (@files[@filenos]) {
	($files[0]{linetype} eq $$href{linetype})
	    or die ("file '$$href{name}', line $$href{lineno}, expected $files[0]{linetype}, ".
		    "found $$href{linetype}\n");
	
	if ($$href{linetype} eq 'marker line') {
	    unless ($options{relax} == 1) {
		($files[0]{mrknum} == $$href{mrknum}) 
		    or die ("file '$$href{name}', line $$href{lineno}, loci number mismatch\n");
		(join (',', @{$files[0]{mrknames}}) eq join (',', @{$$href{mrknames}}))
		    or die ("file '$$href{name}', line $$href{lineno}, loci name mismatch\n");
	    }
	    
	} elsif ($$href{linetype} eq 'header line') {
	    (join (',', sort (keys (%{$files[0]{cmpheaders}}))) eq 
	     join (',', sort (keys (%{$$href{cmpheaders}}))))
		or die ("file '$$href{name}', line $$href{lineno}, header name mismatch\n");
	    print ("headers: ", join (', ', @{$$href{headers}}), "\n");
	    print ("cmpheaders: ", join (', ', keys (%{$$href{cmpheaders}})), "\n");
	    print ("regex: $$href{regex}\n");
	    print ("fmt: $$href{fmt}\n");
	    
	} elsif ($$href{linetype} eq 'data line') {
	    foreach $key (keys (%{$files[0]{cmpheaders}})) {
		(($key eq 'AVG_LR') || ($key eq 'PPL')) and next;
		$va = $files[0]{cmpheaders}{$key};
		($files[0]{flds}[$va] eq $$href{flds}[$va])
		    or die ("file '$$href{name}', line $$href{lineno}, ".
			    "expected $files[0]{flds}[$va] in field '$key'\n");
	    }
	    $va = $files[0]{cmpheaders}{AVG_LR};
	    $files[0]{flds}[$va] *= $$href{flds}[$va];
	}
    }

    if ($files[0]{linetype} eq 'marker line') {
	print ("# $files[0]{mrknum} ", join (' ', @{$files[0]{mrknames}}), "\n");

    } elsif ($files[0]{linetype} eq 'header line') {
	print (join (' ', @{$files[0]{headers}}), "\n");
	    
    } elsif ($files[0]{linetype} eq 'data line') {
	if ($options{mode} eq 'multipoint') {
	    $va = $files[0]{cmpheaders}{AVG_LR};
	    ((($avglr = $files[0]{flds}[$va]) < 0.214) || 
	     (($ppl = ($avglr * $avglr) / (-5.77 + (54 * $avglr) + ($avglr * $avglr))) < 0))
		and $ppl = 0;
	    $va = $files[0]{cmpheaders}{PPL};
	    $files[0]{flds}[$va] = $ppl;
	}
	printf ($files[0]{fmt}, @{$files[0]{flds}});

    } elsif ($files[0]{linetype} eq 'end of file') {
	exit (0);
    }
}


sub get_next_line
{
    my ($f) = @_;
    my $buff;
    my ($mrknum, $mrknames);
    my (@headers, $base, @exts, $ext, $str, @regex, $keepers);
    my (@flds);
    my $va;

    $$f{lineno}++;
    $buff = $$f{fp}->getline;
    if (! defined ($buff)) {
	$$f{linetype} = 'end of file';
	return (1);
    }

    $buff =~ s/^\s*(\S|\S.*\S)\s*$/$1/s;

    if (($mrknum, $mrknames) = ($buff =~ /\#\s+(\d+)((?:\s+\S+)+)/)) {
	$$f{linetype} = 'marker line';
	$$f{expect} = 'header line';
	$mrknames =~ s/^\s*(\S|\S.*\S)\s*$/$1/s;
	$$f{mrknum} = $mrknum;
	@{$$f{mrknames}} = split (/\s+/, $mrknames);

    } elsif ($buff =~ /AVG_?LR/i) {
	$$f{linetype} = 'header line';
	$$f{expect} = 'data line';

	# If the regex field is non-empty, we've seen a header line
	# before, and we save ourselves a lot of work by returning now.
	# This assumes that header line (and consequently data line) 
	# formats won't change in mid-file.

	($$f{regex} ne '')
	    and return (1);
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
		# field in the data line. Build the regex to match the same number
		# of reals.

		@exts = split (/,/, $ext);
		$str = '\(' . '[\-\d\.]+' . '[,\s]+[\-\d\.]+' x (scalar (@exts)-1) . '\)';

	    } else {
		# Otherwise it's just a simple field.

		$base = $headers[$va];
		$str = '[\-\d\.]+';
	    }
	    if ($base =~ /(D\d\d|Dprime|Theta|Pos|Markers)/i) {
		# D-Prime, Theta, or Position fields we need to keep, so build
		# the regex such that those fields will be captured.

		push (@regex, '(' . $str . ')');
		$$f{regex} = ($$f{regex}) ? join ('\s+', $$f{regex}, @regex) :
		    join ('\s+', @regex) ;
		@regex = ();
		$$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' %s' : '%s';
		$$f{cmpheaders}{$headers[$va]} = $keepers;
		$keepers++;
		
	    } elsif ($base =~ /PPL/i) {
		# Same as with D-Prime, etc., except we force the field name 
	    
		push (@regex, '(' . $str . ')');
		$$f{regex} = ($$f{regex}) ? join ('\s+', $$f{regex}, @regex) :
		    join ('\s+', @regex) ;
		@regex = ();
		$$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' %.4f' : '%.4f';
		$$f{cmpheaders}{PPL} = $keepers;
		$keepers++;

	    } elsif ($base =~ /AVG_?LR/i) {
		# Same with Avg_LR field, except we force the header name in
		# the cmpheaders hash to 'AVG_LR' so we can look for that 
		# specific string later. Also, if the mode is multipoint and
		# there is exactly one parenthesized extention to the field
		# name, force the regex to accomodate the 'avg_LR(count)'
		# format: recognize but don't capture the parethesized column.

		if (($options{mode} eq 'multipoint') && (scalar (@exts) == 1)) {
		    push (@regex, '([\-\d\.]+)\(\d+\)');
		    $$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' %.4f(0)' : '%.4f(0)';
		} else {
		    push (@regex, '(' . $str . ')');
		    $$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' %.4f' : '%.4f';
		}
		$$f{regex} = ($$f{regex}) ? join ('\s+', $$f{regex}, @regex) :
		    join ('\s+', @regex) ;
		@regex = ();
		$$f{cmpheaders}{AVG_LR} = $keepers;
		$keepers++;
		
	    } elsif (($base =~ /pen_?vector/i) && ($options{mode} eq 'multipoint')) {
		# An unfortunate hack: multipoint output only emits one header
		# for the three columns of the penetrance vector.

		map { push (@regex, $str); } (0, 1, 2);
		$$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' 0 0 0' : '0 0 0';

	    } elsif (($base =~ /markerlist/i) && ($options{mode} eq 'multipoint')) {
		# Another unfortunate hack: multipoint output emits a plain header
		# for a parenthesized field that contains a comma-seperated list
		# of marker numbers. Force the regex accordingly.

		push (@regex, '(\([\d\,]+\))');
		$$f{regex} = ($$f{regex}) ? join ('\s+', $$f{regex}, @regex) :
		    join ('\s+', @regex) ;
		@regex = ();
		$$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' %s' : '%s';

	    } else {
		# Chaffe fields. Build the regex to recognize them, but not capture.

		push (@regex, $str);
		$$f{fmt} = ($$f{fmt}) ? $$f{fmt} . ' 0' : '0';
	    }
	    push (@{$$f{headers}}, $headers[$va]);
	    $va++;
	}
	$$f{fmt} .= "\n";

    } elsif ($buff =~ /^[\-\d\.]+/) {
	$$f{linetype} = 'data line';
	($$f{regex})
	    or die ("file '$$f{name}', line $$f{lineno}, no previous header line\n");
	(@flds = ($buff =~ /^$$f{regex}/))
	    or die ("file '$$f{name}', line $$f{lineno}, malformed data line\n");
	@{$$f{flds}} = @flds;
    }

    return (1);
}


=head1 USAGE

seq_update_ld.pl [-m mode] [--relax] avgfile1 avgfile2 [avgfile3 ...]

Performs a sequential update of the average likelihood ratios (avgLR)
across multilple files.

=head1 COMMAND LINE ARGUMENTS

=item -m <twopoint|multipoint>

=over 3

Sets the mode. The default is 'twopoint'.

=back

=item --relax

=over 3

Specifies that the loci names in the marker lines need not be the same
in all input files.

=back

=cut
