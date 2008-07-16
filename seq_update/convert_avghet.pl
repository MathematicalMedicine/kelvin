#!perl -w
use strict;

my $line;
my $lineno = 1;
my $version;
my $current = "0.35.0";
my ($maj, $min, $patch);
my $force_chr = '';
my @flds;

# These will be filled in by the appropriate header parsing routine
my $header = '';
my $data_regex = '';
my @fld_idxs = ();

my %header_parsers = ( "unknown" => \&pre_0_35_0,
		       "0.35.0" => \&ver_0_35_0 );

my $usage = "usage: $0 [ -c <chr_num> ] <filename>\n";

if (scalar (@ARGV)) {
    if (($force_chr) = ($ARGV[0] =~ /-c(\d+)?/)) {
	shift (@ARGV);
	($force_chr || ($force_chr = shift (@ARGV))) or die ($usage);
	($force_chr =~ /^\d+$/) or die ($usage);
	$force_chr .= ' ';
    } else {
	$force_chr = '';
    }
}

(defined ($line = <>)) or die ("$0: empty input file\n");

# Version information line, only appears in 0.35.0 and later
if (($version) = ($line =~ /\#\s+Version\s+(\S+)/)) {
    (($maj, $min, $patch) = ($version =~ /(\d+)\.(\d+)\.(\d+)/))
	or die ("$0: badly formed version number '$version' at line $lineno\n");
    $version = sprintf ("%d.%d.%d", $maj, $min, $patch);
    $lineno++;
    (defined ($line = <>))
	or die ("$0: input file ends unexpectedly at line $lineno\n");
} else {
    $version = 'unknown';
}
($version eq $current) 
    and die ("input file is already in most recent format ($current)\n");
exists ($header_parsers{$version})
    or die ("$0: don't know how to parse headers for version '$version'\n");
(($version eq 'unknown') && (length ($force_chr) <= 0))
    and die ("$0: -c <chr_num> is required for this file\n");

print ("# Version $current (converted from $version)\n");

while (1) {

    if ($line =~ /^\#/) {
	# 2-Point marker info line
	print ($line);

    } elsif ($line !~ /^[\-\+\d\.\s\(\,\)e]+$/) {
	# Column header line
	(length ($header) <= 0)
	    and &{$header_parsers{$version}} ($line);
	print ($header);

    } else {
	# By elimination, a data line
	(length ($data_regex) > 0)
	    or die ("$0: data line found before any header lines at line $lineno\n");
	$line =~ s/\s+/ /g;
	$line =~ s/, /,/g;

	(@flds) = ($line =~ /$data_regex/)
	    or die ("$0: badly formed data line at line $lineno\n");
	print ($force_chr, join (' ', @flds[@fld_idxs]), "\n");
    }
    
    (defined ($line = <>)) or last;
    $lineno++;
}


sub pre_0_35_0 {
    my ($line) = @_;
    my @names;
    my @dprimes = ();
    my %mapped = ();
    my $str;
    my $va = 0;

    # Trim leading and trailing whitespace
    $line =~ s/^\s*(\S|\S.*\S)\s*$/$1/;
    # Remove spaces following commas
    $line =~ s/,\s+/,/;
    # And split on whitespace
    @names = split (/\s+/, $line);

    while ($str = shift (@names)) {
	$str = uc ($str);
	if ($str =~ /^D\d\d/) {
	    push (@dprimes, $str);
	    $mapped{$str} = $va;
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "DPRIME") {
	    push (@dprimes, "D11");
	    (exists ($mapped{"D11"}))
		and die ("$0: multiple 'DPrime' columns in header at line $lineno\n");
	    $mapped{D11} = $va;
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "THETA(M,F)") {
	    $mapped{"Theta(M,F)"} = $va;
	    $data_regex .= '(\([\d\.]+,[\d\.]+\)) ';
	    
	} elsif ($str eq "COUNT") {
	    $data_regex .= '(\d+) ';

	} elsif ($str =~ /^(AVG_?LR|BR)$/) {
	    $mapped{"BayesRatio"} = $va;
	    $data_regex .= '([\-\d\.]+(?:[eE][\+\-]\d+)?) '

	} elsif ($str eq "AVGLR(COUNT)") {
	    $mapped{"BayesRatio"} = $va;
	    $data_regex .= '([\-\d\.]+(?:[eE][\+\-]\d+)?)\(\d+\) '

	} elsif (($str eq "MAX_HLOD") || ($str eq "MOD")) {
	    $mapped{"MOD"} = $va;
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "R2") {
	    $mapped{"R2"} = $va;
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "ALPHA") {
	    $mapped{"Alpha"} = $va;
	    $data_regex .= '([\d\.]+) ';

	} elsif ($str eq "DGF") {
	    $mapped{"DGF"} = $va;
	    $data_regex .= '([\d\.]+) ';

	} elsif ($str eq "MF") {
	    $mapped{"MF"} = $va;
	    $data_regex .= '([\d\.]+) ';

	} elsif ($str eq "PEN_DD") {
	    (($names[0] =~ /PEN_DD/i) && ($names[1] =~ /PEN_DD/i))
		or die ("$0: one or more penetrance headers missing at line $lineno\n");
	    splice (@names, 0, 2);
	    $mapped{"PenetranceVector"} = $va;
	    $data_regex .= '([\d\.]+ [\d\.]+ [\d\.]+) ';

	} elsif ($str eq "PEN_VECTOR") {
	    $mapped{"PenetranceVector"} = $va;
	    $data_regex .= '([\d\.]+ [\d\.]+ [\d\.]+) ';

	} elsif ($str =~ /^POS(ITION)?/) {
	    $mapped{"Position"} = $va;
	    $data_regex .= '([\d\.]+) ';
	    
	} elsif ($str eq "PPL") {
	    $mapped{"PPL"} = $va;
	    $data_regex .= '([\d\.]+) ';

	} elsif ($str eq "MARKERLIST") {
	    $mapped{"MarkerList"} = $va;
	    $data_regex .= '(\([\d\,]+\)) ';

	} else {
	    die ("$0: unknown header '$str' at line $lineno\n");
	}
	$va++;
    }
    $header = 'Chr ';
    foreach $str (@dprimes) {
	$header .= "$str ";
	push (@fld_idxs, $mapped{$str});
    }
    foreach $str ('Theta(M,F)', qw/Position PPL BayesRatio MOD R2 Alpha
		  DGF MF PenetranceVector MarkerList/) {
	defined ($mapped{$str}) or next;
	$header .= "$str ";
	push (@fld_idxs, $mapped{$str});
    }
    
    chop ($data_regex, $header);
    $header .= "\n";
    return;
}


# A stub that should never be called (at least, not until 0.36 is released)
sub ver_0_35_0 {
    ;
}
