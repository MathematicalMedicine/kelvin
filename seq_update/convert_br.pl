#!/usr/bin/perl -w
use strict;

#
# convert_br - Convert Bayes Ratio (or avghet.out) files generated by previous
#              versions of Kelvin to the most recent format.
#
# John Burian - john.burian@nationwidechildrens.org
#
# Copyright 2008, The Research Institute at Nationwide Children's Hospital
# All rights reserved. Permission is granted to use this software for
# non-profit educational purposes only.
#

my $line;
my $lineno = 1;
my $version;
my $current = "0.35.0";
my ($maj, $min, $patch);
my @flds;
my @arr;
my $va;

my $force_chr = '';
my $markerlist_len = 0;
my $penvector_type = '';
my $liability_classes = 1;

# These globals are set in the unkown version header parser, and used
# to determine the penetrance vector type, the number of liability 
# classes and the number of markers
my $markerlist_col = -1;
my $penvector_col = -1;

# These will be filled in by the appropriate header parsing routine
my $header = '';
my $data_regex = '';
my $data_fmt = '';
my @fld_idxs = ();

my %header_parsers = ( "unknown" => \&pre_0_35_0,
		       "0.35.0" => \&ver_0_35_0 );

my $usage = "usage: $0 [ -c <chr_num> ] <filename>\n";

if (scalar (@ARGV)) {
    if (($force_chr) = ($ARGV[0] =~ /-c(\d+)?/)) {
        shift (@ARGV);
        ($force_chr || ($force_chr = shift (@ARGV))) or die ($usage);
        ($force_chr =~ /^\d+$/) or die ($usage);
    } elsif ($ARGV[0] =~ /(\?|help)/) {
        system ("/usr/bin/pod2text $0");
        exit (0);
    } else {
        $force_chr = '';
    }
}
(scalar (@ARGV) == 1) or die ($usage);
open (FH, $ARGV[0]) or die ("open '$ARGV[0]' failed, $!|n");
(defined ($line = <FH>)) or die ("$0: empty input file\n");

# Version information line, only appears in 0.35.0 and later
if (($version) = ($line =~ /\#\s+Version\s+(\S+)/)) {
    (($maj, $min, $patch) = ($version =~ /(\d+)\.(\d+)\.(\d+)/))
	or die ("$0: badly formed version number '$version' at line $lineno\n");
    $version = sprintf ("%d.%d.%d", $maj, $min, $patch);
    $lineno++;
    (defined ($line = <FH>))
	or die ("$0: input file ends unexpectedly at line $lineno\n");
} else {
    $version = 'unknown';
}
(versionnum ($version) >= versionnum ($current)) 
    and die ("input file is already in most recent format ($current)\n");
exists ($header_parsers{$version})
    or die ("$0: don't know how to parse headers for version '$version'\n");
(($version eq 'unknown') && (length ($force_chr) <= 0))
    and die ("$0: -c <chr_num> is required for this file\n");

# We have to build a general $data_regex and try to match one data line,
# then count fields in PenetranceVector and MarkerList to figure out the
# penetrance vector type, number of liability classes and how the converted
# header line should look like.
if ($version eq 'unknown') {

    # 2-Point marker info line
    if ($line =~ /^\#/) {
	(defined ($line = <FH>))
	    or die ("$0: input file ends unexpectedly at line $lineno\n");
	$lineno++;
    }
    &{$header_parsers{$version}} ($line);

    (defined ($line = <FH>))
	or die ("$0: input file ends unexpectedly at line $lineno\n");
    $lineno++;
    $line =~ s/\s+/ /g;
    $line =~ s/, +/,/g;
    (@flds) = ($line =~ /^\s*$data_regex\s*$/)
	or die ("$0: badly formed data line at line $lineno\n");
    ($markerlist_col != -1)
	and $markerlist_len = scalar (@arr = split (/,/, $flds[$markerlist_col]));
    $va = scalar (@arr = split (/ /, $flds[$penvector_col]));
    if (($va == 3) || ($va == 6) || ($va == 9)) {
	$penvector_type = "DT";
	$liability_classes = $va / 3;
    } elsif (($va == 4) || ($va == 8)) {
	$penvector_type = "QT_Chi";
	$liability_classes = $va / 4;
    } elsif (($va == 7) || ($va == 14)) {
	$penvector_type = "QT_T";
	$liability_classes = $va / 7;
    } else {
	die ("$0: can't determine penetrance vector type at line $lineno\n");
    }

    # Reset $header to an empty string, which will force the first header line
    # to be re-parsed, informed by $penvector_type, etc.
    $header = '';
    # And rewind the file to the beginning. We ignore the possibility of a 
    # version line, because we wouldn't be here if there had been one.
    seek (FH, 0, 0);
    $line = <FH>;
    $lineno = 1;
}

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
	$line =~ s/, +/,/g;

	(@flds) = ($line =~ /^\s*$data_regex\s*$/)
	    or die ("$0: badly formed data line at line $lineno\n");
	printf ($data_fmt, @flds[@fld_idxs]);
    }
    
    (defined ($line = <FH>)) or last;
    $lineno++;
}


sub pre_0_35_0 {
    my ($line) = @_;
    my @names;
    my @dprimes = ();
    my %mapped = ();
    my $str;
    my $va = 0;
    my $lc;

    $header = 'Chr ';
    $data_fmt = "$force_chr ";
    $data_regex = '';
    @fld_idxs = ();
    # Trim leading and trailing whitespace
    $line =~ s/^\s*(\S|\S.*\S)\s*$/$1/;
    # Remove spaces following commas
    $line =~ s/, +/,/;
    # And split on whitespace
    @names = split (/\s+/, $line);

    while ($str = shift (@names)) {
	$str = uc ($str);
	if ($str =~ /^D\d\d/) {
	    push (@dprimes, $str);
	    $mapped{$str} = [ $va ];
	    $data_regex .= '([\-\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "DPRIME") {
	    push (@dprimes, "D11");
	    (exists ($mapped{"D11"}))
		and die ("$0: multiple 'DPrime' columns in header at line $lineno\n");
	    $mapped{D11} = [ $va ];
	    $data_regex .= '([\-\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "THETA(M,F)") {
	    $mapped{"Theta(M,F)"} = [ $va ];
	    $data_regex .= '(\([\d\.]+,[\d\.]+\)) ';
	    $data_fmt .= '%s ';
	    
	} elsif ($str eq "COUNT") {
	    $data_regex .= '(\d+) ';

	} elsif ($str =~ /^(AVG_?LR|BR)$/) {
	    $mapped{"BayesRatio"} = [ $va ];
	    $data_regex .= '([\-\d\.]+(?:[eE][\+\-]\d+)?) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "AVGLR(COUNT)") {
	    $mapped{"BayesRatio"} = [ $va ];
	    $data_regex .= '([\-\d\.]+(?:[eE][\+\-]\d+)?)\(\d+\) ';
	    $data_fmt .= '%s ';

	} elsif (($str eq "MAX_HLOD") || ($str eq "MOD")) {
	    $mapped{"MOD"} =  [$va ];
	    $data_regex .= '([\-\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "R2") {
	    $mapped{"R2"} = [ $va ];
	    $data_regex .= '([\-\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "ALPHA") {
	    $mapped{"Alpha"} = [ $va ];
	    $data_regex .= '([\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "DGF") {
	    $mapped{"DGF"} = [ $va ];
	    $data_regex .= '([\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif ($str eq "MF") {
	    $mapped{"MF"} = [ $va ];
	    $data_regex .= '([\d\.]+) ';
	    $data_fmt .= '%s ';

	} elsif (($str eq "PEN_VECTOR") || ($str eq "PEN_DD")) {
	    if ($str eq "PEN_DD") {
		(($names[0] =~ /PEN_DD/i) && ($names[1] =~ /PEN_DD/i))
		    or die ("$0: one or more penetrance headers missing at line $lineno\n");
		splice (@names, 0, 2);
	    }
	    $mapped{"PenetranceVector"} = [ ];
	    if ($penvector_type eq "") {
		push (@{$mapped{"PenetranceVector"}}, $va);
		$penvector_col = $va;
		$data_regex .= '([\-\d\.]+(?: [\-\d\.]+)+) ';
		$data_fmt .= '%s ';
	    } elsif ($penvector_type eq "DT") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    push (@{$mapped{"PenetranceVector"}}, $va, $va+1, $va+2);
		    $data_regex .= '([\-\d\.]+) ' x 3;
		    $data_fmt .= '(%s,%s,%s) ';
		    $va += 3;
		}
		$va--;
	    } elsif ($penvector_type eq "QT_Chi") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    push (@{$mapped{"PenetranceVector"}}, $va, $va+1, $va+2, $va+3);
		    $data_regex .= '([\-\d\.]+) ' x 4;
		    $data_fmt .= '(%s,%s,%s,%s) ';
		    $va += 4;
		}
		$va--;
	    } elsif ($penvector_type eq "QT_T") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    push (@{$mapped{"PenetranceVector"}}, $va, $va+1, $va+2, $va+3, $va+4,
			  $va+5, $va+6);
		    $data_regex .= '([\-\d\.]+) ' x 7;
		    $data_fmt .= '(%s,%s,%s,%s,%s,%s,%s) ';
		    $va += 7;
		}
		$va--;
	    }
	    
	} elsif ($str =~ /^POS(ITION)?/) {
	    $mapped{"Position"} = [ $va ];
	    $data_regex .= '([\d\.]+) ';
	    $data_fmt .= '%s ';
	    
	} elsif ($str eq "PPL") {
	    $mapped{"PPL"} = [ $va ];
	    $data_regex .= '([\d\.]+) ';
	    $data_fmt .= '%s ';
	    
	} elsif ($str eq "MARKERLIST") {
	    $mapped{"MarkerList"} = [ $va ];
	    if ($markerlist_len == 0) {
		$markerlist_col = $va;
	    } 
	    $data_regex .= '(\([\d\,]+\)) ';
	    $data_fmt .= '%s ';
	    
	} else {
	    die ("$0: unknown header '$str' at line $lineno\n");
	}
	$va++;
    }
    
    foreach $str (@dprimes) {
	$header .= "$str ";
	push (@fld_idxs, @{$mapped{$str}});
    }
    foreach $str ('Theta(M,F)', qw/Position PPL BayesRatio MOD R2 Alpha
		  DGF MF PenetranceVector MarkerList/) {
	defined ($mapped{$str}) or next;
	if ($str eq "PenetranceVector") {
	    if ($penvector_type eq "") {
		$header .= "$str ";
	    } elsif ($penvector_type eq "DT") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    $header .= "LC${lc}PV(DD,Dd,dd) ";
		}
	    } elsif ($penvector_type eq "QT_Chi") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    $header .= "LC${lc}PV(DDMean,DdMean,ddMean,Thresh) ";
		}
	    } elsif ($penvector_type eq "QT_T") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    $header .= "LC${lc}PV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD,Thresh) ";
		}
	    }
	} elsif ($str eq "MarkerList") {
	    $header .= "MarkerList(". join (',', (0 .. ($markerlist_len - 1))) . ") ";
	} else {
	    $header .= "$str ";
	}
	push (@fld_idxs, @{$mapped{$str}});
    }
    
    chop ($data_regex, $data_fmt, $header);
    $header .= "\n";
    $data_fmt .= "\n";
    return;
}


# A stub that should never be called (at least, not until 0.36 is released)
sub ver_0_35_0
{
    ;
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

convert_br.pl [-c <chrnum>] infile

Converts Bayes ratio (br.out) or average heterogeneity (avghet.out) files
generated by older version of Kelvin the most recent br.out format.

=head1 COMMAND LINE ARGUMENTS

=item -c <chrnum>

=over 3

Specifiy the chromosome number. Required if the input file does not contain
a chromosome number column, and optional otherwise. If specified, and the input
file already contains a chromosome column, the chromosome number in the output
will be forced to <chrnum>.

=back

=cut
