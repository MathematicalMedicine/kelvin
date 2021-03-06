#!/usr/bin/env perl
use strict;
use warnings;

#
# convert_br - Convert Bayes Ratio (or avghet.out) files generated by previous
#              versions of Kelvin to the most recent format.
#
# John Burian - john.burian@nationwidechildrens.org
#
# Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
#

my $line;
my $lineno = 1;
my $version;
my $current = "0.38.0";
my ($maj, $min, $patch);
my $br_file = undef;
my @flds;
my @arr;
my $va;

my $force_chr = '';
# These are all used when a map file is necessary to insert a position 
# column in the converted file.
my $map_file = '';
my $map_fmt = undef;
my ($mrk_idx, $pos_idx);
my %marker_pos = ();

# How many markers in a multipoint file
my $markerlist_len = 0;
# Penetrance vector format: DT, QT w/ std dist, or QT w/ Chi**2 dist
my $penvector_type = '';
my $liability_classes = 1;

# These globals are set in the unkown version header parser, and used
# to determine the penetrance vector type, the number of liability 
# classes and the number of markers
my $markerlist_col = -1;
my $penvector_col = -1;


# This will be filled in by the appropriate marker info parsing routine
my %marker_flds;
# This will be filled in by the appropriate header parsing routine
my %header_info;

my $twopoint = 0;

my %header_parsers = ( "unknown" => \&pre_0_35_0,
		       "0.35.0" => \&ver_0_35_0,
		       "0.36.0" => \&ver_0_35_0,
		       "0.36.1" => \&ver_0_36_1,
		       "0.37.0" => \&ver_0_36_1,
		       "0.37.1" => \&ver_0_36_1,
		       "0.37.2" => \&ver_0_36_1,
		       "0.37.3" => \&ver_0_36_1,
		       "0.37.4" => \&ver_0_36_1,
		       "0.37.5" => \&ver_0_36_1,
		       "0.37.6" => \&ver_0_36_1,
		       "0.37.7" => \&ver_0_36_1,
		       "0.38.0" => \&ver_0_38_0 );

my $usage = "usage: $0 [-c <chr_num>] [-m <map_file>] <filename>\n";

while (scalar (@ARGV)) {
    if ($ARGV[0] =~ /^-c(\d+)?/) {
        shift (@ARGV);
	$force_chr = (defined ($1)) ? $1 : shift (@ARGV);
        (defined ($force_chr) && ($force_chr =~ /^\d+$/)) or die ($usage);
    } elsif ($ARGV[0] =~ /^-m(\S+)?/) {
        shift (@ARGV);
	$map_file = (defined ($1)) ? $1 : shift (@ARGV);
        (defined ($map_file)) or die ($usage);
    } elsif ($ARGV[0] =~ /(\?|help)/) {
        system ("/usr/bin/pod2text $0");
        exit (0);
    } else {
	defined ($br_file) and die ($usage);
	$br_file = shift (@ARGV);
    }
}

if (($map_file) && ($map_file ne 'nomap')) {
    open (FH, $map_file) or die ("$0: map file '$map_file': $!\n");
    while (! defined ($map_fmt)) {
	(defined ($line = <FH>)) or die ("$0: unknown map format in '$map_file'\n");
	if ($line =~ /\s*CHR(OMOSOME)?\s+(HALDANE|KOSAMBI|POSITION)\s+(NAME|MARKER)/i) {
	    $map_fmt = '\d+\s+([\-\d\.]+(?:e[\-\+]\d+)?)\s+(\S+)';
	    ($mrk_idx, $pos_idx) = (1, 0);
	} elsif ($line =~ /\s*CHR(OMOSOME)?\s+(NAME|MARKER)\s+(HALDANE|KOSAMBI|POSITION)/i) {
	    $map_fmt = '\d+\s+(\S+)\s+([\-\d\.]+(?:e[\-\+]\d+)?)';
	    ($mrk_idx, $pos_idx) = (0, 1);
	} elsif ($line =~ /\d+\s+([\-\d\.]+(?:e[\-\+]\d+)?)\s+(\S+)/) {
	    $marker_pos{$2} = $1;
	    $map_fmt = '\d+\s+([\-\d\.]+(?:e[\-\+]\d+)?)\s+(\S+)';
	    ($mrk_idx, $pos_idx) = (1, 0);
	} elsif ($line =~ /\d+\s+(\S+)\s+([\-\d\.]+(?:e[\-\+]\d+)?)/) {
	    $marker_pos{$1} = $2;
	    $map_fmt = '\d+\s+(\S+)\s+([\-\d\.]+(?:e[\-\+]\d+)?)';
	    ($mrk_idx, $pos_idx) = (0, 1);
	}
    }
    while ($line = <FH>) {
	unless (@arr = ($line =~ /$map_fmt/)) {
	    chomp ($line);
	    die ("$0: badly formatted line in '$map_file': '$line'\n");
	}
	$marker_pos{$arr[$mrk_idx]} = $arr[$pos_idx];
    }
    close (FH);
}

open (FH, $br_file) or die ("open '$br_file' failed, $!\n");
(defined ($line = <FH>)) or die ("$0: empty input file\n");

# Version information line, only appears in 0.35.0 and later
# All files with no Version line are treated as 'unknown'
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
if (versionnum ($version) >= versionnum ($current)) {
    # and die ("input file is already in most recent format ($current)\n");
    # A cheap hack: don't die, just copy out files that are already current
    print ("# Version V$version\n", $line, <FH>);
    exit (0);
}
exists ($header_parsers{$version})
    or die ("$0: don't know how to parse headers for version '$version'\n");

if ($version eq 'unknown') {
    # Unversioned files, two-point and multipoint both, 
    (length ($force_chr) <= 0) and die ("$0: -c <chr_num> is required for this file\n");
    
    # We have to build a general $data_regex and try to match one data line,
    # then count fields in PenetranceVector and MarkerList to figure out the
    # penetrance vector type, number of liability classes and how the converted
    # header line should look like.

    # A 2-point marker info line. Unversioned 2-pt files require a map
    if ($line =~ /^\#/) {
	($map_file) or die ("$0: can't convert without a map file (use -m <map_file>)\n");
	$lineno++;
	$twopoint = 1;
	(defined ($line = <FH>))
	    or die ("$0: input file ends unexpectedly at line $lineno\n");
    }
    &{$header_parsers{$version}} ($line, \%header_info);

    (defined ($line = <FH>))
	or die ("$0: input file ends unexpectedly at line $lineno\n");
    $lineno++;
    $line =~ s/\s+/ /g;
    $line =~ s/, +/,/g;

    (@flds) = ($line =~ /^\s*$header_info{data_regex}\s*$/)
	or die ("$0: badly formed data line at line $lineno\n");
    ($markerlist_col != -1)
	and $markerlist_len = scalar (@arr = split (/,/, $flds[$markerlist_col]));
    $va = scalar (@arr = split (/ /, $flds[$penvector_col]));

    # Imprinting wasn't supported before 0.36, so we don't need to consider how
    # that changes the penetrance vector length. Also note that for all
    # unversioned files, QT penetrance vectors included a Threshold column, 
    # even if the run was straight QT and not CT (this is true even in some
    # versioned files, but at least those are labeled).

    if (($va == 3) || ($va == 6) || ($va == 9)) {
	# Dichotomous trait: three elements (DD Dd dd) per liability class
	$penvector_type = "DT";
	$liability_classes = $va / 3;
    } elsif (($va == 4) || ($va == 8)) {
	# Chi-Squared: Three degrees-of-freedom (DD Dd dd) and a threshold 
	# per liability class
	$penvector_type = "QT_Chi";
	$liability_classes = $va / 4;
    } elsif (($va == 7) || ($va == 14)) {
	# Normal or T: Three means (DD Dd dd) and three standard deviations (DD Dd dd)
	# and a threshold per liability class.
	$penvector_type = "QT_T";
	$liability_classes = $va / 7;
    } else {
	die ("$0: can't determine penetrance vector type at line $lineno\n");
    }

    # Rewind the file to the beginning. We ignore the possibility of a 
    # version line, because we wouldn't be here if there had been one.
    seek (FH, 0, 0);
    $line = <FH>;
    $lineno = 1;

} elsif (versionnum ($version) < versionnum ("0.36.1")) {
    # Versioned 2-point files less than V0.36.1 have no position column, so
    # we require a map. 

    if ($line =~ /^\#/) {
	($map_file) or die ("$0: can't convert without a map file (use -m <map_file>)\n");
	$twopoint = 1;
    }
}

print ("# Version V$current (converted from ". (($version eq "unknown") ? "" : "V"). "$version)\n");
while (1) {

    if ($line =~ /^\#/) {
	# 2-Point marker info line
	$twopoint = 1;
	%marker_flds = ();
	if (($version eq 'unknown') || (versionnum ($version) < versionnum ("0.38.0"))) {
	    pre_0_38_0_marker ($line, \%marker_flds);
	} else {
	    ver_0_38_0_marker ($line, \%marker_flds);
	}

    } elsif ($line !~ /^[\-\+\d\.\s\(\,\)e]+$/) {
	# Column header line
	&{$header_parsers{$version}} ($line, \%header_info);

    } else {
	# By elimination, a data line
	(length ($header_info{data_regex}) > 0)
	    or die ("$0: data line found before any header lines at line $lineno\n");
	$line =~ s/\s+/ /g;
	$line =~ s/, +/,/g;
	(@flds) = ($line =~ /^\s*$header_info{data_regex}\s*$/)
	    or die ("$0: badly formed data line at line $lineno\n");

	(exists ($header_info{chr_col}) && ! exists ($marker_flds{Chr}))
	    and $marker_flds{Chr} = $flds[$header_info{chr_col}];
	(exists ($header_info{pos_col}) && !
	 (exists ($marker_flds{Position}) || exists ($marker_flds{Position1}) ||
	  exists ($marker_flds{AvgPosition})))
	    and $marker_flds{Position} = $flds[$header_info{pos_col}];

	if (length ($header_info{header}) > 0) {
	    # Guard against Perl thinking 'build_marker' is a file handle
	    ($twopoint) and print ("". build_marker (\%marker_flds));
	    print ($header_info{header});
	    $header_info{header} = '';
	}
	printf ($header_info{data_fmt}, @flds[@{$header_info{fld_idxs}}]);
    }
    
    (defined ($line = <FH>)) or last;
    $lineno++;
}


sub pre_0_38_0_marker
{
    my ($line, $href) = @_;
    my ($seq, $name1, $name2);

    (($seq, $name1, $name2) = ($line =~ /\#\s+(\d+)\s+(\S+)\s+(\S+)/))
	or die ("$0: badly formatted marker info at line $lineno\n");
    $$href{Seq} = $seq;
    if (! $map_file) {
	$$href{Trait} = $name1;
	$$href{Marker} = $name2;

    } elsif ($map_file eq 'nomap') {
	$$href{Trait} = $name1;
	$$href{Marker} = $name2;
	$$href{Position} = $seq;

    } elsif (! exists ($marker_pos{$name1})) {
	$$href{Trait} = $name1;
	(exists ($marker_pos{$name2}))
	    or die ("$0: can't find position based on marker info at line $lineno\n");
	$$href{Marker} = $name2;
	$$href{Position} = $marker_pos{$name2};

    } elsif (exists ($marker_pos{$name2}))  {
	$$href{Marker1} = $name1;
	$$href{Position1} = $marker_pos{$name1};
	$$href{Marker2} = $name2;
	$$href{Position2} = $marker_pos{$name2};
    } else {
	die ("$0: can't find position based on marker info at line $lineno\n");
    }
    ($force_chr) and $$href{Chr} = $force_chr;

    return (1);
}


# This should never be called, but we implement it now while the details are
# fresh in mind, in case we need it later.
sub ver_0_38_0_marker
{
    my ($line, $href) = @_;
    my ($fld, $val);
    my %legal = (Seq => '', Chr => '', Trait => '', Marker => '', Position => '',
		 AvgPosition => '', MalePosition => '', FemalePosition => '',
		 Marker1 => '', Position1 => '', Marker2 => '', Position2 => '',
		 Physical => '');

    $line =~ s/^\#\s*//;
    while ($line !~ /^\s*$/) {
	($fld, $val, $line) = split (/\s+/, $line, 3);
	(($fld =~ s/:$//) && exists ($legal{$fld}))
	    or die ("$0: unknown field '$fld' in marker info at line $lineno\n");
	$$href{$fld} = $val;
    }
    return (1);
}


sub build_marker
{
    my ($href) = @_;
    my $line = "#";

    map { $line .= " $_: $$href{$_}" } qw(Seq Chr);
    if (exists $$href{Position}) {
	map { $line .= " $_: $$href{$_}" } qw(Trait Marker Position);
    } elsif (exists ($$href{AvgPosition})) {
	map {
	    $line .= " $_: $$href{$_}";
	} qw(Trait Marker AvgPosition MalePosition FemalePosition);
    } elsif (exists ($$href{Position1})) {
	map { $line .= " $_: $$href{$_}" } qw(Marker1 Position1 Marker2 Position2);
    }
    (exists ($$href{Physical})) and $line .= " Physical: $$href{Physical}";
    $line .= "\n";

    return ($line);
}


sub pre_0_35_0
{
    my ($line, $href) = @_;
    my @names;
    my @dprimes = ();
    my %mapped = ();
    my $str;
    my $va = 0;
    my $lc;
    my $header = '';
    my $data_regex = '';
    my $data_fmt = '';
    my $fld_idxs = [];

    %$href = ();
    if ($twopoint) {
	$header = $data_fmt = '';
    } else {
	$header = 'Chr ';
	$data_fmt = "$force_chr ";
    }
    $data_regex = '';
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
	    $mapped{$str}{idx} = [ $va ];
	    $mapped{$str}{fmt} = '%.2f ';
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "DPRIME") {
	    push (@dprimes, "D11");
	    (exists ($mapped{"D11"}))
		and die ("$0: multiple 'DPrime' columns in header at line $lineno\n");
	    $mapped{D11}{idx} = [ $va ];
	    $mapped{D11}{fmt} = '%.2f ';
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "THETA(M,F)") {
	    $mapped{"Theta(M,F)"}{idx} = [ $va ];
	    $mapped{"Theta(M,F)"}{fmt} = '%s ';
	    $data_regex .= '(\([\d\.]+,[\d\.]+\)) ';
	    
	} elsif ($str eq "COUNT") {
	    $data_regex .= '(\d+) ';

	} elsif ($str =~ /^(AVG_?LR|BR|BAYESRATIO)$/) {
	    $mapped{"BayesRatio"}{idx} = [ $va ];
	    $mapped{"BayesRatio"}{fmt} = '%e ';
	    $data_regex .= '([\-\d\.]+(?:[eE][\+\-]\d+)?) ';

	} elsif ($str eq "AVGLR(COUNT)") {
	    $mapped{"BayesRatio"}{idx} = [ $va ];
	    $mapped{"BayesRatio"}{fmt} = '%s ';
	    $data_regex .= '([\-\d\.]+(?:[eE][\+\-]\d+)?)\(\d+\) ';

	} elsif (($str eq "MAX_HLOD") || ($str eq "MOD")) {
	    $mapped{"MOD"}{idx} =  [$va ];
	    $mapped{"MOD"}{fmt} = '%.4f ';
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "R2") {
	    $mapped{"R2"}{idx} = [ $va ];
	    $mapped{"R2"}{fmt} = '%.4f ';
	    $data_regex .= '([\-\d\.]+) ';

	} elsif ($str eq "ALPHA") {
	    $mapped{"Alpha"}{idx} = [ $va ];
	    $mapped{"Alpha"}{fmt} = '%.2f ';
	    $data_regex .= '([\d\.]+) ';

	} elsif ($str eq "DGF") {
	    $mapped{"DGF"}{idx} = [ $va ];
	    $mapped{"DGF"}{fmt} = '%.4f ';
	    $data_regex .= '([\d\.]+) ';

	} elsif ($str eq "MF") {
	    $mapped{"MF"}{idx} = [ $va ];
	    $mapped{"MF"}{fmt} = '%.4f ';
	    $data_regex .= '([\d\.]+) ';

	} elsif (($str eq "PEN_VECTOR") || ($str eq "PEN_DD")) {
	    if ($str eq "PEN_DD") {
		(($names[0] =~ /PEN_DD/i) && ($names[1] =~ /PEN_DD/i))
		    or die ("$0: one or more penetrance headers missing at line $lineno\n");
		splice (@names, 0, 2);
	    }
	    $mapped{"PenetranceVector"}{idx} = [ ];
	    $mapped{"PenetranceVector"}{fmt} = '';
	    if ($penvector_type eq "") {
		push (@{$mapped{"PenetranceVector"}{idx}}, $va);
		$mapped{"PenetranceVector"}{fmt} .= '%s ';
		$penvector_col = $va;
		$data_regex .= '([\-\d\.]+(?: [\-\d\.]+)+) ';
	    } elsif ($penvector_type eq "DT") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    push (@{$mapped{"PenetranceVector"}{idx}}, $va, $va+1, $va+2);
		    $mapped{"PenetranceVector"}{fmt} .= '(%s,%s,%s) ';
		    $data_regex .= '([\-\d\.]+) ' x 3;
		    $va += 3;
		}
		$va--;
	    } elsif ($penvector_type eq "QT_Chi") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    push (@{$mapped{"PenetranceVector"}{idx}}, $va, $va+1, $va+2, $va+3);
		    $mapped{"PenetranceVector"}{fmt} .= '(%s,%s,%s,%s) ';
		    $data_regex .= '([\-\d\.]+) ' x 4;
		    $va += 4;
		}
		$va--;
	    } elsif ($penvector_type eq "QT_T") {
		for ($lc = 0; $lc < $liability_classes; $lc++) {
		    push (@{$mapped{"PenetranceVector"}{idx}}, $va, $va+1, $va+2, $va+3,
			  $va+4, $va+5, $va+6);
		    $mapped{"PenetranceVector"}{fmt} .= '(%s,%s,%s,%s,%s,%s,%s) ';
		    $data_regex .= '([\-\d\.]+) ' x 7;
		    $va += 7;
		}
		$va--;
	    }
	    
	} elsif ($str =~ /^POS(ITION)?/) {
	    $mapped{"Position"}{idx} = [ $va ];
	    $mapped{"Position"}{fmt} = '%.4f ';
	    $data_regex .= '([\d\.]+) ';
	    
	} elsif ($str eq "PPL") {
	    $mapped{"PPL"}{idx} = [ $va ];
	    $mapped{"PPL"}{fmt} = '%s ';
	    $data_regex .= '([\d\.]+) ';
	    
	} elsif ($str eq "MARKERLIST") {
	    $mapped{"MarkerList"}{idx} = [ $va ];
	    $mapped{"MarkerList"}{fmt} = '%s ';
	    if ($markerlist_len == 0) {
		$markerlist_col = $va;
	    } 
	    $data_regex .= '(\([\d\,]+\)) ';
	    
	} else {
	    die ("$0: unknown header '$str' at line $lineno\n");
	}
	$va++;
    }

    foreach $str (@dprimes) {
	$header .= "$str ";
	push (@$fld_idxs, @{$mapped{$str}{idx}});
	$data_fmt .= $mapped{$str}{fmt};
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
	push (@$fld_idxs, @{$mapped{$str}{idx}});
	$data_fmt .= $mapped{$str}{fmt};
    }
    
    chop ($data_regex, $data_fmt, $header);
    $header .= "\n";
    $data_fmt .= "\n";

    $$href{header} = $header;
    $$href{data_regex} = $data_regex;
    $$href{data_fmt} = $data_fmt;
    $$href{fld_idxs} = $fld_idxs;
    return;
}


# V0.35.0 through V0.36.0 : Remove the Chr column from two-point files
sub ver_0_35_0
{
    my ($line, $href) = @_;

    %$href = ();
    if (! $twopoint) {
	# Multipoint. No changes needed, really.
	$$href{header} = $line;
	$$href{data_regex} = '(\S.*\S)';
	$$href{data_fmt} = "%s\n";
	$$href{fld_idxs} = [0];
    } else {
	# Otherwise, twopoint. We need to strip the Chr column from the data lines,
	# but still capture the Chr value in the regex so we can insert the Chr value
	# in the marker info line.
	($$href{header} = $line) =~ s/Chr\s+(\S.*)/$1/;
	$$href{data_regex} = '(\d+)\s+(\S.*)';
	$$href{data_fmt} = "%s\n";
	$$href{fld_idxs} = [1];
	$$href{chr_col} = 0;
    }
    return;
}


# V0.36.1 and V0.37.x: Remove the Chr and Position columns from two-point files
sub ver_0_36_1
{
    my ($line, $href) = @_;

    %$href = ();
    if (! $twopoint) {
	# Multipoint. No changes needed, really.
	$$href{header} = $line;
	$$href{data_regex} = '(\S.*\S)';
	$$href{data_fmt} = "%s\n";
	$$href{fld_idxs} = [0];
    } else {
	# Otherwise, twopoint. We need to strip the Chr and Position columns
	# from the data lines, but still capture those values in the regex so
	# we can insert them in the marker info line.
	($$href{header} = $line) =~ s/Chr\s+Position\s+(\S.*)/$1/;
	$$href{data_regex} = '(\d+)\s+([\d\.]+)\s+(\S.*)';
	$$href{data_fmt} = "%s\n";
	$$href{fld_idxs} = [2];
	$$href{chr_col} = 0;
	$$href{pos_col} = 1;
    }
    return;
}


# A stub that should never be called
sub ver_0_38_0
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


sub lookup_pos
{
    my ($line) = @_;
    my ($num, $name1, $name2);

    unless (($num, $name1, $name2) = ($line =~ /\#\s+(\d+)\s+(\S+)\s+(\S+)/)) {
	chomp ($line);
	die ("$0: badly formatted marker info at line $lineno\n");
    }
    ($map_file eq 'nomap') and return ($num);
    exists ($marker_pos{$name1}) and return ($marker_pos{$name1});
    exists ($marker_pos{$name2}) and return ($marker_pos{$name2});
    die ("$0: can't find position based on marker info at line $lineno\n");
}


=head1 USAGE

convert_br.pl [-c <chr_num>] [-m <map_file>] infile

Converts Bayes ratio (br.out) or average heterogeneity (avghet.out) files
generated by older version of Kelvin the most recent br.out format.

=head1 COMMAND LINE ARGUMENTS

=item -c <chr_num>

=over 3

Specifies the chromosome number. Required if the input file does not contain
a chromosome number column, and optional otherwise. If specified, and the input
file already contains a chromosome column, the chromosome number in the output
will be forced to <chrnum>.

=back

=item -m <map_file>

=over 3

Specifies a map file that lists centiMorgan positions of markers. Required if
the input file contains two-point Bayes ratios and does not contain a marker
position column, and ignored otherwise. If map_file is specified as the literal
string 'nomap', the marker positions in the output will be filled with the
sequential marker number.

=back

=cut
