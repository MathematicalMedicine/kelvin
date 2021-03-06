#!/usr/bin/env perl
# Copyright (C) 2009, 2022 Mathematical Medicine LLC
# 
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
use strict;
use warnings;

my $line;
my %state = (points => '',
	     trait => '',
	     distrib => '',
	     grid => '',
	     poly => '',
	     classes => 1,
	     imprinting => 0,
	     librium => '',
	     sexspecific => 0,
	     marker2marker => 0
	     );
my @phenocodes = ();
my @traitpositions = ();
my @thetas = ();
my @malethetas = ();
my @femalethetas = ();
my @dprimes = ();
my @genefreqs = ();
my @allelefreqs = ();
my @alphas = ();
my @params = ();
my @thresholds = ();
my %pens = ( DD => [], Dd => [], dD => [], dd => [] );
my %constraints = (pensimple => [], penclass => [], paramsimple => [], paramclass => []);
my $brfile = undef;
my $pplfile = undef;
my @arr;

my @files = ();
my @modelparam = ();
my @trait = ();
my @type = ();
my @logs = ();

while ($line = <>) {
    $line =~ s/\#.*$//;
    ($line =~ /^\s*$/) and next;
    $line =~ s/^\s*(\S|\S.*\S)\s*$/$1/;
    $line =~ s/\s+/ /g;

    if ($line =~ /^SA (\d+)$/) {
	push (@type, "Multipoint $1\n");
	$state{points} = 'multi';

    } elsif ($line =~ /^SS (\d+)$/) {
	push (@type, "Multipoint $1\n");
	push (@type, "SexSpecific\n");
	$state{points} = 'multi';
	$state{sexspecific} = 1;

    } elsif ($line =~ /^SS$/) {
	push (@type, "SexSpecific\n");
	$state{sexspecific} = 1;

    } elsif ($line =~ /^TP$/) {
	push (@type, "# Twopoint is the default\n");
	$state{points} = 'two';

    } elsif ($line =~ /^IMP$/) {
	push (@type, "Imprinting\n");
	$state{imprinting} = 1;

    } elsif ($line =~ /^MM$/) {
	push (@type, "MarkerToMarker all\n");
	$state{marker2marker} = 1;

    } elsif ($line =~ /^AM$/) {
	push (@type, "MarkerToMarker adjacent\n");
	$state{marker2marker} = 1;

    } elsif ($line =~ /^LC (\d+)$/) {
	if ($1 > 1) {
	    $state{classes} = $1;
	    push (@type, "LiabilityClasses $1\n");
	}

    } elsif ($line =~ /^PE *(\d+)?/) {
	print ("# Polynomial Evaluation is the default\n");
	(defined ($1) && ($1 > 1)) and push (@type, "PolynomialScale $1\n");
	$state{poly} = 'on';

    } elsif ($line =~ /^DK$/) {
	print ("# Dynamic sampling is the default\n");
	$state{grid} = 'dynamic';


    } elsif ($line =~ /^DT$/) {
	push (@trait, "# Dichotomous Trait is the default\n");
	$state{trait} = 'DT';

    } elsif ($line =~ /^QT\s+normal\s+(\S+)\s+(\S+)/) {
	$state{trait} = 'QT';
	$state{distrib} = "normal $1, $2";

    } elsif ($line =~ /^QT\s+T\s+(\S+)\s+(\S+)\s+(\S+)/) {
	$state{trait} = 'QT';
	$state{distrib} = "normal $2, $3";

    } elsif ($line =~ /^QT\s+chisq\s+(\S+)\s+(\S+)/) {
	$state{trait} = 'QT';
	$state{distrib} = 'chisq';

    } elsif ($line =~ /^TT (\S.*)$/) {
	push (@thresholds, parserange ($1));

    } elsif ($line =~ /^XC$/) {
	push (@type, "SexLinked\n");

    } elsif ($line =~ /^UP (\d+)$/) {
	push (@type, "# Unknown person code is fixed at 0\n");

    } elsif ($line =~ /^AS (\S+) (\S+) (\S+)$/) {
	push (@phenocodes, $1, $2, $3);

    } elsif ($line =~ /^DA (\d+)$/) {
	push (@type, "DiseaseAlleles $1\n");


    } elsif ($line =~ /^LD (\S.*)$/) {
	push (@dprimes, parserange ($1));
	$state{points} = 'two';
	$state{librium} = 'disequi';

    } elsif ($line =~ /^Th (\S.*)$/) {
	push (@thetas, parserange ($1));

    } elsif ($line =~ /^Tm (\S.*)$/) {
	push (@malethetas, parserange ($1));
	$state{sexspecific} = 1;

    } elsif ($line =~ /^Tf (\S.*)$/) {
	push (@femalethetas, parserange ($1));
	$state{sexspecific} = 1;

    } elsif ($line =~ /^TL (\S.*)$/) {
	push (@traitpositions, parserange ($1));

    } elsif ($line =~ /^TM$/) {
	push (@traitpositions, "Marker");

    } elsif ($line =~ /^GF (\S.*)$/) {
	push (@genefreqs, parserange ($1));

    } elsif ($line =~ /^AF (\S.*)$/) {
	push (@allelefreqs, parserange ($1));

    } elsif ($line =~ /^AL (\S.*)$/) {
	push (@alphas, parserange ($1));

    } elsif ($line =~ s/^\s*(dd)\s*(>=|>|!=|==)\s*(dd)(\s*;)?//i) {
	push (@{$constraints{pensimple}}, "$1 $2 $3");
	while ($line =~ s/^\s*(dd)\s*(>=|>|!=|==)\s*(dd)(\s*;)?//i) {
	    $constraints{pensimple}[-1] .= ", $1 $2 $3";
	}

    } elsif ($line =~ s/^\s*(dd)\s*(\d+)\s*(>=|>|!=|==)\s*(dd)\s*(\d+)(\s*;)?//i) {
	push (@{$constraints{penclass}}, "$1 $2 $3 $4 $5");
	while ($line =~ s/^\s*(dd)\s*(\d+)\s*(>=|>|!=|==)\s*(dd)\s*(\d+)(\s*;)?//i) {
	    $constraints{penclass}[-1] .= ", $1 $2 $3 $4 $5";
	}

    } elsif ($line =~ s/^\s*P1\s*(dd)\s*(>=|>|!=|==)\s*P1\s*(dd)(\s*;)?//i) {
	push (@{$constraints{paramsimple}}, "$1 $2 $3");
	while ($line =~ s/^\s*P1\s*(dd)\s*(>=|>|!=|==)\s*P1\s*(dd)(\s*;)?//i) {
	    $constraints{paramsimple}[-1] .= ", $1 $2 $3";
	}

    } elsif ($line =~ s/^\s*P1\s*(dd)\s*(\d+)\s*(>=|>|!=|==)\s*P1\s*(dd)\s*(\d+)(\s*;)?//i) {
	push (@{$constraints{paramclass}}, "$1 $2 $3 $4 $5");
	while ($line =~ s/^\s*P1\s*(dd)\s*(\d+)\s*(>=|>|!=|==)\s*P1\s*(dd)\s*(\d+)(\s*;)?//i) {
	    $constraints{paramclass}[-1] .= ", $1 $2 $3 $4 $5";
	}

    } elsif ($line =~ /^(dd) (\S.*)$/i) {
	push (@{$pens{$1}}, parserange ($2));

    } elsif ($line =~ s/^P1 (\S.*)/$1/i) {
	# A little hokey-pokey because P1 can be specified differently
	$line =~ s/[\s;]+/;/g;
	push (@params, parserange ($line));

    } elsif ($line =~ /^PD (\S+)/) {
	push (@files, "PedigreeFile $1\n");

    } elsif ($line =~ /^DF (\S+)/) {
	push (@files, "LocusFile $1\n");

    } elsif ($line =~ /^MK (\S+)/) {
	push (@files, "FrequencyFile $1\n");

    } elsif ($line =~ /^MP (\S+)/) {
	push (@files, "MapFile $1\n");

    } elsif ($line =~ /^CF (\S+)/) {
	push (@files, "CountFile $1\n");

    } elsif ($line =~ /^MD (\S+)/) {
	push (@files, "MODFile $1\n");

    } elsif ($line =~ /^MX (\S+)/) {
	push (@files, "ExtraMODs\n");

    } elsif ($line =~ /^DIR (\S+)/) {
	push (@files, "NIDetailFile $1\n");

    } elsif ($line =~ /^HE (\S+)/) {
	$brfile = $1;
	      
    } elsif ($line =~ /^PF (\S+)/) {
	$pplfile = $1;
	      
    } elsif ($line =~ /^LOG/i) {
	push (@logs, $line. "\n");

    } elsif ($line =~ /^(T_MIN|T_MAX)/) {
	;

    } else {
	print ("not handled: '$line'\n");
    }
}

if (scalar (@{$pens{dD}}) && ! $state{imprinting}) {
    $state{imprinting} = 1;
    push (@type, "Imprinting\n");
}

@type = sort (@type);
@trait = sort (@trait);
@modelparam = sort (@modelparam);
@files = sort (@files);

($state{trait} eq '') and $state{trait} = 'DT';
($state{points} eq '') and $state{points} = 'two';
($state{grid} eq '') and $state{grid} = 'fixed';
($state{librium} eq '') and $state{librium} = 'equi';
($state{poly} eq '') and $state{poly} = 'off';


(($state{points} eq 'multi') && (scalar (@traitpositions)))
    and push (@type, "TraitPositions ". join (", ", @traitpositions). "\n");
($state{librium} eq 'disequi')
    and push (@type, "LD\n");
($state{poly} eq 'off') and push (@type, "NonPolynomial\n");

if ($state{grid} eq 'fixed') {
    push (@modelparam, "FixedModels\n");
    if ($state{points} eq 'two') {
	if ($state{sexspecific}) {
	    push (@type, "SexSpecific\n");
	    (scalar (@malethetas))
		and push (@modelparam, "MaleTheta ". join (", ", @malethetas). "\n");
	    (scalar (@femalethetas))
		and push (@modelparam, "FemaleTheta ". join (", ", @femalethetas). "\n");
	} else {
	    (scalar (@thetas)) and push (@modelparam, "Theta ". join (", ", @thetas). "\n");
	}
	(($state{librium} eq 'disequi') && (scalar (@dprimes)))
	    and push (@modelparam, "DPrime ". join (", ", @dprimes). "\n");
    }
    if (! $state{marker2marker}) {

	(scalar (@genefreqs))
	    and push (@modelparam, "DiseaseGeneFrequency ". join (", ", @genefreqs). "\n");
	(scalar (@alphas)) and push (@modelparam, "Alpha ". join (", ", @alphas). "\n");
	
	if ($state{trait} eq 'DT') {
	    foreach (qw(DD Dd dD dd)) {
		(exists ($pens{$_}) && scalar (@{$pens{$_}})) or next;
		push (@modelparam, "Penetrance $_ ". join (", ", @{$pens{$_}}). "\n");
	    }
	    if (scalar (@{$constraints{pensimple}})) {
		map { push (@modelparam, "Constrain Penetrance $_\n") } @{$constraints{pensimple}};
	    }
	    if (($state{classes} > 1) && (scalar (@{$constraints{penclass}}))) {
		map { push (@modelparam, "Constrain Penetrance $_\n") } @{$constraints{penclass}};
	    }
	    
	} elsif ($state{trait} eq 'QT') {
	    (scalar (@thresholds)) and $state{trait} = 'QTT';
	    push (@trait, "$state{trait} $state{distrib}\n");
	    (scalar (@thresholds)) and push (@trait, "Threshold ". join (", ", @thresholds). "\n");
	    if ($state{distrib} =~ /normal/i) {
		foreach (qw(DD Dd dD dd)) {
		    (exists ($pens{$_}) && scalar (@{$pens{$_}})) or next;
		    push (@modelparam, "Mean $_ ". join (", ", @{$pens{$_}}). "\n");
		}
		(scalar (@params))
		    and push (@modelparam, "StandardDev ". join (", ", @params). "\n");
		if (scalar (@{$constraints{pensimple}})) {
		    map { push (@modelparam, "Constrain Mean $_\n") } @{$constraints{pensimple}};
		}
		if (($state{classes} > 1) && (scalar (@{$constraints{penclass}}))) {
		    map { push (@modelparam, "Constrain Mean $_\n") } @{$constraints{penclass}};
		}
		if (scalar (@{$constraints{paramsimple}})) {
		    map { push (@modelparam, "Constrain StandardDev $_\n")
			  } @{$constraints{paramsimple}};
		}
		if (($state{classes} > 1) && (scalar (@{$constraints{paramclass}}))) {
		    map { push (@modelparam, "Constrain StandardDev $_\n")
			  } @{$constraints{paramclass}};
		}
		
	    } elsif ($state{distrib} =~ /chisq/) {
		foreach (qw(DD Dd dD dd)) {
		    (exists ($pens{$_}) && scalar (@{$pens{$_}})) or next;
		    push (@modelparam, "DegreesOfFreedom $_ ". join (", ", @{$pens{$_}}). "\n");
		}
		if (scalar (@{$constraints{pensimple}})) {
		    map { push (@modelparam, "Constrain DegreesOfFreedom $_\n")
			  } @{$constraints{pensimple}};
		}
		if (($state{classes} > 1) && (scalar (@{$constraints{penclass}}))) {
		    map { push (@modelparam, "Constrain DegreesOfFreedom $_\n")
			  } @{$constraints{penclass}};
		}
	    }
	}
    }
    (scalar (@allelefreqs))
        and push (@modelparam, "MarkerAlleleFrequency ". join (", ", @allelefreqs). "\n");


} else {
    if ($state{trait} eq 'QT') {
	(scalar (@thresholds)) and $state{trait} = 'QTT';
	push (@trait, "$state{trait} $state{distrib}\n");
	if (scalar (@thresholds)) {
	    @arr = expandrange (@thresholds);
	    push (@trait, "Threshold $arr[0], $arr[-1]\n");
	}
	if ($state{distrib} =~ /chisq/i) {
	    foreach (qw(DD Dd dD dd)) {
		(exists ($pens{$_}) && scalar (@{$pens{$_}})) or next;
		@arr = expandrange (@{$pens{$_}});
		push (@modelparam, "DegreesOfFreedom $_ $arr[0], $arr[-1]\n");
	    }
	}
    }
    ($state{sexspecific}) and push (@type, "SexSpecific\n");

}

if (scalar (@phenocodes) && ! $state{marker2marker}) {
    if ($state{trait} ne 'QT') {
	push (@trait, "PhenoCodes ". join (", ", @phenocodes). "\n");
    } else {
	push (@trait, "PhenoCodes $phenocodes[0]\n");
    }
}


(($state{points} eq 'two') && defined ($pplfile))
    and push (@files, "PPLFile $pplfile\n");

# This is how we really should behave WRT marker-to-marker and BR files
#(defined ($brfile) && ! $state{marker2marker})
#    and push (@files, "BayesRatioFile $brfile\n");

# This is how we're going to behave for now, to make LaVonne's life easier
(defined ($brfile))
    and push (@files, "BayesRatioFile $brfile\n");
($state{marker2marker})
    and push (@files, "ForceBRFile\n");

# push (@files, "MODFile mod.out\n");
# ($state{points} eq 'two')
#     and push (@files, "ExtraMODs\n");


unshift (@type, "# Type\n");
unshift (@trait, "# Trait\n");
unshift (@modelparam, "# Modelparam\n");
unshift (@files, "# Files\n");

print (@type, "\n");
print (@trait, "\n");
print (@modelparam, "\n");
print (@files, "\n");
(scalar (@logs)) and print (@logs);

sub parserange
{
    my ($range) = @_;
    my @arr;

    if ($range =~ /^([\-\d\.]+) ([\-\d\.]+) ([\d\.]+)$/) {
	push (@arr, "$1-$2:$3");

    } elsif ($range =~ /^([\-\d\.]+)$/) {
	push (@arr, $1);

    } elsif ($range =~ /^([\-\d\.]+)\s*;/) {
	@arr = split (/[\s;]+/, $range);

    } else {
	die ("parserange: can't handle '$range'\n");
    }
    return (@arr);
}


sub expandrange
{
    my @range = @_;
    my $el;
    my ($val, $start, $end, $incr);
    my $va;
    my @vals = ();

    foreach $el (@range) {
	if ($el =~ /^[\-\d\.]+$/) {
	    push (@vals, $el);
	} elsif (($start, $end, $incr) = ($el =~ /([\-\d\.]+)-([\-\d\.]+):([\-\d\.]+)/)) {
	    for ($va = 0; ($val = $start + $va * $incr) <= $end; $va++) {
		push (@vals, $val);
	    }
	} else {
	    die ("expandrange: can't handle '$el'\n");
	}
    }
    @vals = sort ({$a <=> $b} @vals);
    return (@vals);
}
