#!/usr/bin/env perl
# Copyright (C) 2010, 2022 Mathematical Medicine LLC
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

#
# KelvinFamily: an object for managing groups of related KelvinIndividuals (families).
#
package KelvinFamily v1.10.0;
our $errstr='';

sub new
{
    my ($class, $aref) = @_;
    my $family = bless ({}, $class);
    my $dataset;
    my $makeped = undef;
    my $phased = undef;
    my $traits;
    my $trait = undef;
    my $affection;
    my $multigeneration = 0;
    my ($aff_founders, $unaff_founders, $aff_kids, $unaff_kids) = (0, 0, 0, 0);
    my ($dad, $mom, $kids) = (undef, undef, []);
    my $href = {};
    my $va;

    @$family{qw/pedid pedtype count founders nonfounders/} = (undef, undef, 0, 0, 0);
    @$family{qw/founderpairs multmarriages individuals/} = (0, 0, []);
    @$family{qw/dad mom kids origfmt loopids/} = (undef, undef, undef, undef, undef);
    
    unless (scalar (@$aref) > 0) {
	$errstr = "array of individuals is empty";
	return (undef);
    }
    $dataset = $$aref[0]{dataset};
    $makeped = $$aref[0]{makeped};
    $phased = $$aref[0]{phased};
    $$family{pedid} = $$aref[0]{pedid};
    $$family{count} = scalar (@$aref);

    $$href{$$aref[0]{indid}} = '';
    for ($va = 1; $va < $$family{count}; $va++) {
	if ($$aref[$va]{dataset} != $dataset) {
	    $errstr = "individuals are based on different datasets";
	    return (undef);
	}
	if ($$aref[$va]{makeped} ne $makeped) {
	    $errstr = "individuals are mixed pre- and post-makeped format";
	    return (undef);
	}
	if ($$aref[$va]{phased} ne $phased) {
	    $errstr = "individuals are mixed phased and non-phased genotypes";
	    return (undef);
	}
	if ($$aref[$va]{pedid} ne $$family{pedid}) {
	    $errstr = "individuals are from different families (pedigree IDs)";
	    return (undef);
	}
	if (exists ($$href{$$aref[$va]{indid}})) {
	    $errstr = "individual $$aref[$va]{indid} appears more than once";
	    return (undef);
	}
	$$href{$$aref[$va]{indid}} = '';
    }
    
    @{$$family{individuals}} = @$aref;
    # Sort individuals such that no one precedes their parents. This is
    # required for makeped(), and handy for find_loops().
    if ($$family{count} > 1) {
        $family->sort_family or return (undef);
    }

    if ($$aref[0]{makeped} eq 'pre') {
	$family->makeped or return (undef);
	$$family{origfmt} = 'pre';
    } else {
	($$family{count} == 1 || $family->verify_links) or return (undef);
	$$family{origfmt} = 'post';
    }
    $family->find_loops;
    if (defined ($$family{loopids}) && $$family{origfmt} eq "post") {
        $errstr = "family $$family{pedid} has loops involving individuals $$family{loopids}";
        return (undef);
    }
    
    (defined ($traits = $dataset->traitOrder))
	and $trait = $$traits[0];
    for ($va = 0; $va < $$family{count}; $va++) {
        # TODO: shouldn't assume default affection status coding (0, 1, 2)
	$affection = (defined ($trait)) ? $$aref[$va]->getTrait ($trait) : 1;
	if ($$aref[$va]{dadid} eq '0' && $$aref[$va]{momid} eq '0') {
	    $$family{founders}++;
	    if ($$aref[$va]{sex} == 1) {
		$dad = $$aref[$va];
	    } else {
		$mom = $$aref[$va];
	    }
	    if ($affection eq '2') {
		$aff_founders++;
	    } else {
		$unaff_founders++;
	    }
	} else {
	    ($$aref[$va]{firstchildid} ne '0')
		and $multigeneration = 1;
	    $$family{nonfounders}++;
	    push (@$kids, $$aref[$va]);
	    if ($affection eq '2') {
		$aff_kids++;
	    } else {
		$unaff_kids++;
	    }
	}
    }

    if ($$family{count} == 1 && $$family{nonfounders} == 0) {
	# TODO: additional verification? require undefined parents? require genotype?
	$$family{pedtype} = 'casecontrol';
	$$family{kids} = [ defined ($dad) ? $dad : $mom ];
	$$family{nonfounders} = 1;
	$$family{founders} = 0;
    } elsif ($$family{count} == 3 && $$family{founders} == 2) {
	if (! ($$dad{genotyped} || $$dad{phenotyped} || $$mom{genotyped} || $$mom{phenotyped})) {
	    $$family{pedtype} = 'casecontrol';
	} elsif ($unaff_founders == 2 && $aff_kids == 1) {
	    $$family{pedtype} = 'trio-strict';
	} else {
	    $$family{pedtype} = 'trio';
	}
	@$family{qw/dad mom kids/} = ($dad, $mom, $kids);
    } elsif ($$family{count} == 4 && $$family{founders} == 2) {
	if ($aff_kids == 2) {
	    if ($unaff_founders == 2) {
		$$family{pedtype} = 'asp-strict';
	    } else {
		$$family{pedtype} = 'asp';
            }
	} else {
	    $$family{pedtype} = 'nuclear';
	}
	@$family{qw/dad mom kids/} = ($dad, $mom, $kids);
    } elsif ($$family{founders} < 2 || $$family{count} < 3) {
	$errstr = "family $$family{pedid} has invalid structure";
	return (undef);
    } else {
	$$family{pedtype} = ($$family{founders} > 2 || $multigeneration) ? 'general' : 'nuclear';
	($$family{pedtype} ne 'general')
	    and @$family{qw/dad mom kids/} = ($dad, $mom, $kids);
    }
    return ($family);
}

sub new_from_count
{
    my ($class, $dataset, $pedid, $traits, $labelmap, $individuals) = @_;
    my $family = bless ({}, $class);
    my ($dad, $mom, $kids) = (undef, undef, []);
    my $ind;
    my $va;

    @$family{qw/pedid pedtype count founders nonfounders/} = ($pedid, undef, 0, 0, 0);
    @$family{qw/founderpairs multmarriages individuals/} = (0, 0, []);
    @$family{qw/dad mom kids origfmt/} = (undef, undef, undef, 'post');

    $$family{count} = scalar (@$individuals);
    if ($$family{count} == 1) {
	@{$$family{individuals}} = (KelvinIndividual->new_from_count ($dataset, $pedid, $traits, $labelmap, $$individuals[0]));
	$$family{pedtype} = 'casecontrol';
	$$family{kids} = [ $$family{individuals}[0] ];
	$$family{nonfounders} = 1;
	$$family{founders} = 0;

    } elsif ($$family{count} >= 3) {
	($dad, $mom, @$kids) = map { KelvinIndividual->new_from_count ($dataset, $pedid, $traits, $labelmap, $_) } @$individuals;
	if ($$dad{sex} == 2 && $$mom{sex} == 1) {
	    $ind = $dad;
	    $dad = $mom;
	    $mom = $ind;
	} elsif ($$dad{sex} != 1 || $$mom{sex} != 2) {
	    # TODO: set errstr and return failure
	    die ("illegal sexes in parents");
	}
	# TODO: properly categorize trios and ASPs
	$$family{pedtype} = 'nuclear';
	$$dad{proband} = 1;
	$$mom{indid} = $$mom{origindid} = 2;
	$$family{founders} = 2;
	$$family{nonfounders} = scalar (@$kids);
	$$dad{firstchildid} = $$mom{firstchildid} = 3;
	for ($va = 0; $va < scalar (@$kids); $va++) {
	    $$kids[$va]{indid} = $$kids[$va]{origindid} = $va + 3;
	    $$kids[$va]{patsibid} = $$kids[$va]{matsibid} = $va + 4;
	    @{$$kids[$va]}{qw/dadid momid/} = (1, 2);
	}
	$$kids[-1]{patsibid} = $$kids[-1]{matsibid} = 0;
	@{$$family{individuals}} = ($dad, $mom, @$kids);
	@$family{qw/dad mom kids/} = ($dad, $mom, $kids);
    } else {
	# TODO: set errstr and return failure
	die ("illegal family structure (2 people is enough for love, but not for a pedigree)");
    }
    return ($family);
}

sub map
{
    my ($self, $newset) = @_;
    my $new = {};
    my %hash;
    my $va;
    
    map {
	$$new{$_} = $$self{$_};
    } qw/pedid pedtype count founders nonfounders founderpairs multmarriages origfmt/;

    $$new{individuals} = [];
    for ($va = 0; $va < $$self{count}; $va++) {
	$$new{individuals}[$va] = $$self{individuals}[$va]->map ($newset);
    }

    # Case/control can be a single individual, or it can be one kid with two dummy parents
    if ($$self{pedtype} eq 'casecontrol' && $$self{count} == 1) {
	$$new{dad} = $$new{mom} = undef;
	$$new{kids} = [ $$new{individuals}[0] ];
    } elsif ($$self{pedtype} ne 'general') {
	$hash{$$self{dad}{indid}} = 'dad';
	$hash{$$self{mom}{indid}} = 'mom';
	for ($va = 0; $va < scalar (@{$$self{kids}}); $va++) {
	    $hash{$$self{kids}[$va]{indid}} = $va;
	}
	$$new{kids} = [];
	for ($va = 0; $va < $$new{count}; $va++) {
	    if ($hash{$$new{individuals}[$va]{indid}} eq 'dad') {
		$$new{dad} = $$new{individuals}[$va]; 
	    } elsif ($hash{$$new{individuals}[$va]{indid}} eq 'mom') {
		$$new{mom} = $$new{individuals}[$va];
	    } else {
		$$new{kids}[$hash{$$new{individuals}[$va]{indid}}] = $$new{individuals}[$va];
	    }
	}
    } else {
	$$new{dad} = $$new{mom} = $$new{kids} = undef;
    }
    return (bless ($new, ref ($self)));
}

sub write
{
    my ($self) = @_;
    my $individual;

    foreach $individual (@{$$self{individuals}}) {
	$individual->write or return (undef);
    }
    return (1);
}

sub dataset
{
    my ($self) = @_;

    return ($$self{individuals}[0]{dataset});
}

# Verify the parental, offspring and sibling links in a post-makeped family.
sub verify_links
{
    my ($self) = @_;
    my %hash;
    my $ind;
    my ($pedid, $indid, $dadid, $momid, $firstchildid, $patsibid, $matsibid, $sex, $origindid);
    my $setlist = [];
    my $origpedid;

    $origpedid = $$self{individuals}[0]{origpedid};
    foreach $ind (@{$$self{individuals}}) {
	($pedid, $indid, $dadid, $momid, $firstchildid, $patsibid, $matsibid, $sex, $origindid) = 
	    @$ind{qw/pedid indid dadid momid firstchildid patsibid matsibid sex origindid/};
	foreach ($indid, $firstchildid, $patsibid, $matsibid) {
	    ($_ ne 0 && ! exists ($hash{$_}))
		and $hash{$_} = {sex => undef, dadid => undef, momid => undef, origindid => undef,
				 prevpatsib => undef, prevmatsib => undef, isfirstchild => undef};
	}
	if (defined ($hash{$indid}{sex})) {
	    $errstr = "pedid $pedid, person $indid appears more than once";
	    return (undef);
	}
	@{$hash{$indid}}{qw/dadid momid origindid sex/} = ($dadid, $momid, $origindid, $sex);

	($firstchildid ne '0') and $hash{$firstchildid}{isfirstchild} = 1;
	($patsibid ne '0') and $hash{$patsibid}{prevpatsib} = $indid;
	($matsibid ne '0') and $hash{$matsibid}{prevmatsib} = $indid;
    }

    foreach $indid (keys (%hash)) {
	# If sex is undefined, then the individual ID was inserted due to being
	# referenced as a paternal or maternal sib, or first child, but no such
	# individual was actually present in the pedigree.
	defined ($hash{$indid}{sex}) || delete ($hash{$indid});
    }

    foreach $ind (@{$$self{individuals}}) {
	($pedid, $indid) = @$ind{qw/pedid indid/};
	if ($$ind{dadid} eq '0') {
	    # Founder. KelvinIndividual guarantees that if dad is undefined, mom is also
	    if ($$ind{firstchildid} eq '0') {
		$errstr = "pedid $pedid, person $indid is a founder with no children";
		return (undef);
	    } elsif (! exists ($hash{$$ind{firstchildid}})) {
		$errstr = "pedid $pedid, person $indid missing child $$ind{firstchildid}";
		return (undef);
	    }
	    add_to_set ($setlist, $indid);
	} else {
	    # Non-founder
	    if (! exists ($hash{$$ind{dadid}})) {
		$errstr = "pedid $pedid, person $indid missing father $$ind{dadid}";
		return (undef);
	    } elsif ($$ind{dadid} eq $indid) {
		$errstr = "pedid $pedid, person $indid is coded as their own father";
		return (undef);
	    } elsif ($hash{$$ind{dadid}}{sex} != 1) {
		$errstr = "pedid $pedid, person $indid father $$ind{dadid} is not coded as male";
		return (undef);
	    }
	    if (! exists ($hash{$$ind{momid}})) {
		$errstr = "pedid $pedid, person $indid missing mother $$ind{momid}";
		return (undef);
	    } elsif ($$ind{momid} eq $indid) {
		$errstr = "pedid $pedid, person $indid is coded as it's own mother";
		return (undef);
	    } elsif ($hash{$$ind{momid}}{sex} != 2) {
		$errstr = "pedid $pedid, person $indid mother $$ind{momid} is not coded as female";
		return (undef);
	    }
	    if (! $hash{$indid}{isfirstchild}) {
		if (! defined ($hash{$indid}{prevpatsib})) {
		    $errstr = "pedid $pedid, person $indid is not coded as either a first child or a paternal sibling to other children";
		    return (undef);
		}
		if (! defined ($hash{$indid}{prevmatsib})) {
		    $errstr = "pedid $pedid, person $indid is not coded as either a first child or a maternal sibling to other children";
		    return (undef);
		}
	    }
	    add_to_set ($setlist, @$ind{qw/indid dadid momid/});
	}
	if ($$ind{firstchildid} eq $indid) {
	    $errstr = "pedid $pedid, person $indid is coded as their own child";
	    return (undef);
	}
	if ($$ind{patsibid} ne '0') {
	    if (! exists ($hash{$$ind{patsibid}})) {
		$errstr = "pedid $pedid, person $indid missing sibling $$ind{patsibid}";
		return (undef);
	    } elsif ($$ind{patsibid} eq $indid) {
		$errstr = "pedid $pedid, person $indid is coded as their own sibling";
		return (undef);
	    } elsif ($$ind{dadid} ne $hash{$$ind{patsibid}}{dadid}) {
		$errstr = "pedid $pedid, person $indid and paternal sibling $$ind{patsibid} have different fathers";
		return (undef);
	    }
	}
	if ($$ind{matsibid} ne '0') {
	    if (! exists ($hash{$$ind{matsibid}})) {
		$errstr = "pedid $pedid, person $indid missing sibling $$ind{matsibid}";
		return (undef);
	    } elsif ($$ind{matsibid} eq $indid) {
		$errstr = "pedid $pedid, person $indid is coded is it's own sibling";
		return (undef);
	    } elsif ($$ind{momid} ne $hash{$$ind{matsibid}}{momid}) {
		$errstr = "pedid $pedid, person $indid and maternal sibling $$ind{matsibid} have different mothers";
		return (undef);
	    }
	}
    }

    # $setlist is a ref to an array of sets of connected individuals. There should
    # only be one set, otherwise the family is not completely connected.
    if (scalar (@$setlist > 1)) {
	$errstr = sprintf ("in family %s (orig %s), individuals(s) %s are disconnected from the rest of the pedigree", $$self{pedid}, $origpedid, join (', ', map { "$_ ($hash{$_}{origindid})" } keys (%{$$setlist[-1]})));
	return (undef);
    }

    return (1);
}

# This is an incomplete implementation of MAKEPED. It doesn't handle loops, and it
# always sets the proband to the first founder. Other than that, it should handle
# general pedigrees, multiple marriages, multiple founder pairs, etc.
sub makeped
{
    my ($self) = @_;
    my %hash;
    my $va;
    my $individual;
    my ($pedid, $indid, $dadid, $momid, $kidid);
    my $setlist = [];

    for ($va = 0; $va < scalar (@{$$self{individuals}}); $va++) {
	$individual = $$self{individuals}[$va];
	if ($$individual{makeped} ne 'pre') {
	    $errstr = "ped $pedid mixes pre- and post-MAKEPED formats";
	    return (undef);
	}
	($pedid, $indid, $dadid, $momid) = @$individual{qw/pedid indid dadid momid/};

	if (exists ($hash{$indid})) {
	    $errstr = "pedid $pedid, person $indid appears more than once";
	    return (undef);
	}
	$hash{$indid} = { newid => $va+1, kids => [0], sex => $$individual{sex} };
	if ($dadid eq '0') { 
	    @{$hash{$indid}}{qw/founder patkididx matkididx/} = (1, undef, undef);
	    add_to_set ($setlist, $indid);
	} else {
	    if (! exists ($hash{$dadid})) {
		$errstr = "ped $pedid, person $indid missing father $dadid";
		return (undef);
	    }
	    if (! exists ($hash{$momid})) {
		$errstr = "ped $pedid, person $indid missing mother $momid";
		return (undef);
	    }
	    if ($hash{$dadid}{sex} != 1) {
		$errstr = "ped $pedid, person $indid father $dadid is not coded as male";
		return (undef);
	    }
	    if ($hash{$momid}{sex} != 2) {
		$errstr = "ped $pedid, person $indid mother $momid is not coded as female";
		return (undef);
	    }
	    $hash{$indid}{founder} = 0;
	    $hash{$indid}{patkididx} = scalar (@{$hash{$dadid}{kids}});
	    $hash{$indid}{matkididx} = scalar (@{$hash{$momid}{kids}});
	    # Build kid list with new, rather then current, IDs, to save lookups later
	    # splice (@{$hash{$dadid}{kids}}, -1, 0, $hash{$indid}{newid});
	    # splice (@{$hash{$momid}{kids}}, -1, 0, $hash{$indid}{newid});
	    splice (@{$hash{$dadid}{kids}}, -1, 0, $indid);
	    splice (@{$hash{$momid}{kids}}, -1, 0, $indid);
	    add_to_set ($setlist, $indid, $dadid, $momid);
	}
    }
    # $setlist is a ref to an array of sets of connected individuals. There should
    # only be one set, otherwise the family is not completely connected.
    if (scalar (@$setlist > 1)) {
	$errstr = "in family $$self{pedid}, individual(s) ". join (', ', keys %{$$setlist[-1]}). " disconnected from the rest of the pedigree";
	return (undef);
    }

    for ($va = scalar (@{$$self{individuals}}) - 1; $va >= 0 ; $va--) {
	$individual = $$self{individuals}[$va];
	($indid, $dadid, $momid) = @$individual{qw/indid dadid momid/};
	
	if (! $hash{$indid}{founder}) {
	    $$individual{patsibid} = $hash{$dadid}{kids}[$hash{$indid}{patkididx}];
	    $$individual{matsibid} = $hash{$momid}{kids}[$hash{$indid}{matkididx}];
	    # $$individual{dadid} = $hash{$dadid}{newid};
	    # $$individual{momid} = $hash{$momid}{newid};
	}
	(scalar (@{$hash{$indid}{kids}}))
	    and $$individual{firstchildid} = $hash{$indid}{kids}[0];
	# $$individual{indid} = $hash{$indid}{newid};
    }
    $$self{individuals}[0]{proband} = 1;
    return (1);
}

sub add_to_set
{
    my ($list, @ids) = @_;
    my $merge = 0;
    my $id;
    my ($va, $vb);
    

    for ($va = 0; $va < scalar (@$list); $va++) {
	foreach $id (@ids) {
	    if (exists $$list[$va]{$id})  {
		map { $$list[$va]{$_} = '' } @ids;
		$merge = 1;
	    }
	}
    }
    if (! $merge) {
	push (@$list, {});
	map { $$list[-1]{$_} = '' } @ids;
    } else {
	$va = 0;
	while ($va < scalar (@$list)) {
	    $vb = $va + 1;
	    while ($vb < scalar (@$list)) {
		$merge = 0;
		map { exists ($$list[$va]{$_}) and $merge = 1; } (keys (%{$$list[$vb]}));
		if ($merge) {
		    map { $$list[$va]{$_} = ''; } (keys (%{$$list[$vb]}));
		    splice (@$list, $vb, 1);
		} else {
		    $vb++;
		}
	    }
	    $va++;
	}
    }
    @$list = sort {scalar (keys (%$b)) <=> scalar (keys (%$a))} (@$list);
    return (1);
}

sub sort_family
{
    my ($self) = @_;
    my %parents;
    my %descendants;
    my %famhash;
    my $individual;
    my $indid;
    my (@arr, $aref);
    my $va;

    # Sort a family, such that all founders come first, then non-founders. Inside 
    # each subset, we sort by number of descendents, from greatest to least. This
    # guarantees that no individual appears before their parent.

    # First, a list of 'uncategorized' individuals. As we determine the total number
    # of descendents for each individual, they'll be removed from this list. Everybody
    # starts with 1 'decsendent', representing themselves. We also make a hash of
    # indiviuals, indexed by individual ID, to make it easier to find moms and dads.

    foreach $individual (@{$$self{individuals}}) {
	push (@arr, $$individual{indid});
	if (exists ($famhash{$$individual{indid}})) {
	    $errstr = "pedid $$individual{pedid}, person $$individual{indid} appears more than once";
	    return (undef);
	}
	$famhash{$$individual{indid}} = $individual;
	$descendants{$$individual{indid}} = 1;
    }
    
    # Now iterate until the uncategorized list is empty.
    
    while (scalar (@arr)) {
	
	# Create a hash, the keys of which is the union of the mom and dad IDs for
	# all the individuals in our 'uncategorized' list.

	%parents = ();
	for ($va = 0; $va < scalar (@arr); $va++) {
	    $parents{$famhash{$arr[$va]}{dadid}} = '';
	    $parents{$famhash{$arr[$va]}{momid}} = '';
	}
	# Step through the uncategorized list. If an individual ID does not exist
	# as a key in the hash of mom and dad IDs, then that individual has no 
	# children left in the uncategorized list. We add the individual's
	# descendent count to the descendent counts of both the individual's 
	# parents, then remove the individual from the uncategorized list.

	$va = 0;
	while ($va < scalar (@arr)) {
	    if (exists ($parents{$arr[$va]})) {
		$va++;
	    } else {
		$descendants{$famhash{$arr[$va]}{dadid}} += $descendants{$arr[$va]};
		$descendants{$famhash{$arr[$va]}{momid}} += $descendants{$arr[$va]};
		splice (@arr, $va, 1);
	    }
	}
    }

    # Last, run the family through a sort algorithm that sorts by founder/non-founder,
    # then by number of descendents.

    @{$$self{individuals}} = sort ({by_founder_desc ($a, $b, \%descendants)} values (%famhash));
    return (1);
}

sub by_founder_desc
{
    my ($a, $b, $desc) = @_;
    my $ret;

    if ($$a{dadid} eq '0' && $$a{momid} eq '0') {
	if ($$b{dadid} eq '0' && $$b{momid} eq '0') {
	    (($ret = $$desc{$$b{indid}} <=> $$desc{$$a{indid}}) != 0)
		and return ($ret);
	    (($ret = $$a{sex} <=> $$b{sex}) != 0)
		and return ($ret);
	    return ($$a{indid} cmp $$b{indid});
	}
	return (-1);
    } elsif ($$b{dadid} eq '0' && $$b{momid} eq '0') {
	return (1);
    } elsif (($ret = $$desc{$$b{indid}} <=> $$desc{$$a{indid}}) != 0) {
	return ($ret);
    } elsif (($ret = $$desc{$$b{dadid}} <=> $$desc{$$a{dadid}}) != 0) {
	return ($ret);
    } elsif (($ret = $$desc{$$b{momid}} <=> $$desc{$$a{momid}}) != 0)  {	
	return ($ret);
    } else {
	return ($$a{indid} cmp $$b{indid});
    }
}

sub find_loops
{
    # The logic here is stolen from the kelvin binary. Each individual and
    # mating pair is a node (vertex) in a graph; edges run between the mating
    # pair nodes and the individual nodes for the mates in, or the offspring of,
    # the pair. So a trio (Dad ID 1, Mom ID 2, Offspring ID 3) has four nodes,
    # labeled "1", "2" and "3" (the individuals) and "1+2" (the mating pair).
    # Edges run from "1" to "1+2", "2" to "1+2" and "3" to "1+2". Once the
    # graph is built, we repeatedly traverse the graph, looking for (and
    # removing) nodes with zero or one edge. If we reach a point where no more
    # nodes can be removed, then the individuals remaining are involved in a 
    # loop of some sort (either consanguinous or marriage).

    my ($self) = @_;
    my %nodes;
    my $individual;
    my $key;
    my $lastcount;

    # Nodes are identified by the keys of the hash. Edges are elements in
    # the array refs that are the values of the hash. The edge is recorded
    # in the array ref for both of the nodes that the edge connects.
    foreach $individual (@{$$self{individuals}}) {
        $nodes{$$individual{indid}} = [];
        $$individual{dadid} eq "0" and next;
        $key = $$individual{dadid} . "+" . $$individual{momid};
        if (! exists ($nodes{$key})) {
            push (@{$nodes{$$individual{dadid}}}, $key);
            push (@{$nodes{$$individual{momid}}}, $key);
            $nodes{$key} =  [ $$individual{dadid}, $$individual{momid} ];
        }
        push (@{$nodes{$$individual{indid}}}, $key);
        push (@{$nodes{$key}}, $$individual{indid});
    }

    while ($lastcount = scalar (keys (%nodes))) {
        foreach $key (keys (%nodes)) {
            if (scalar (@{$nodes{$key}}) == 1) {
                # Node has one edge. Delete node and edge.
                delete_edge (\%nodes, $nodes{$key}[0], $key);
                delete ($nodes{$key});
            } elsif (scalar (@{$nodes{$key}}) == 0) {
                # Node has no edges. Delete the node.
                delete ($nodes{$key});
            }
        }
        # True if no nodes were deleted this time through.
        ($lastcount == scalar (keys (%nodes))) and last;
    }
    if ($lastcount > 0) {
        # One or more loops exist. We only want to report the individuals,
        # not the mating pairs.
        foreach $key (keys (%nodes)) {
            $key =~ /\+/ and delete ($nodes{$key});
        }
        $$self{loopids} = join (", ", keys (%nodes));
    }
    return (1);
}

sub delete_edge
{
    my ($nodes, $node, $edge) = @_;
    my $va = 0;

    while ($va < scalar (@{$$nodes{$node}})) {
        if ($$nodes{$node}[$va] eq $edge) {
            splice (@{$$nodes{$node}}, $va, 1);
            return (1);
        }
        $va++;
    }
    return (undef);
}

sub pedid
{
    my ($self) = @_;

    return ($$self{pedid});
}

sub pedtype
{
    my ($self) = @_;

    return ($$self{pedtype});
}

sub origfmt
{
    my ($self) = @_;

    return ($$self{origfmt});
}

sub dad
{
    my ($self) = @_;

    return ($$self{dad});
}

sub mom
{
    my ($self) = @_;

    return ($$self{mom});
}

sub loopids
{
    my ($self) = @_;

    return ($$self{loopids});
}

sub children
{
    my ($self) = @_;
    my $aref = [];
    
    @$aref = (@{$$self{kids}});
    return ($aref);
}

sub individuals
{
    my ($self) = @_;
    my $aref = [];

    @$aref = (@{$$self{individuals}});
    return ($aref);
}


#
# KelvinIndividual: an object for managing a individuals from pedigree files. Note
# that many methods of this object assume intimate knowledge of the internals of the
# KelvinDataset object. Maybe not such great programming practice, but cuts down on
# little utility methods that copy data back and forth.
#
package KelvinIndividual;
our $errstr='';
our $VERSION=1.9;

sub new
{
    my ($class, $dataset, $line) = @_;
    my $trait;
    my $marker;
    my @arr;
    my $traitcol;
    my $markercol;
    my $idx;
    my $colcount;
    my $ind = { pedid => undef, indid => undef, dadid => undef, momid => undef,
		firstchildid => undef, patsibid => undef, matsibid => undef,
		origpedid => undef, origindid => undef, sex => undef,
                proband => undef, traits => [], markers => [], genotyped => 0,
                phenotyped => undef, heterozygous => 0, phased => undef,
                dataset => $dataset, makeped => undef };

    $colcount = scalar (@{$$dataset{markerorder}}) * 2 + scalar (@{$$dataset{traitorder}});
    (defined ($$dataset{undefpheno})) and $$ind{phenotyped} = 0;

    # Pipes indicate phased genotypes. Rearrange the whitespace around the pipes so
    # splitting on whitespace still works like it does with unphased genotypes.
    ($$ind{phased} = ($line =~ s/\s*\|\s*/| /g)) or $$ind{phased} = 0;
    
    # Cut off the original pedigree and person IDs, split on whitespace

    if ($line =~ s/\s*Ped:\s*(\S+)\s+Per:\s*(\S+)\s*$//) {
	# These extra fields are present in genuine post-MAKEPED(tm)-brand pedigree files
	# (accept no substitutes).
	@$ind{qw/origpedid origindid/} = ($1, $2);
	(@$ind{qw/pedid indid dadid momid firstchildid patsibid matsibid sex proband/}, @arr) = 
	    split (' ', $line);
	$$ind{makeped} = 'post';
	if (scalar (@arr) != $colcount) {
	    $errstr = "$$dataset{pedigreefile}, line $$dataset{pedlineno}: too ". ((scalar (@arr) < $colcount) ? "few" : "many"). " columns in pedigree file";
	    return (undef);
	}
	
    } else {
	@arr = split (' ', $line);
	if (scalar (@arr) == $colcount + 9) {
	    (@$ind{qw/pedid indid dadid momid firstchildid patsibid matsibid sex proband/}) =
		splice (@arr, 0, 9);
	    @$ind{qw/origpedid origindid/} = @$ind{qw/pedid indid/};
	    $$ind{makeped} = 'post';
	} elsif (scalar (@arr) == $colcount + 5) {
	    (@$ind{qw/pedid indid dadid momid sex/}) = splice (@arr, 0, 5);
	    @$ind{qw/origpedid origindid/} = @$ind{qw/pedid indid/};
	    @$ind{qw/firstchildid patsibid matsibid proband/} = (0, 0, 0, 0);
	    $$ind{makeped} = 'pre';
	} else {
	    $errstr = "$$dataset{pedigreefile}, line $$dataset{pedlineno}: unexpected number of columns, can't guess format";
	    return (undef);
	}
    }

    if ($$ind{pedid} eq '0') {
	$errstr = "$$dataset{pedigreefile}, $line $$dataset{pedlineno}: illegal pedigree ID '0'";
	return undef;
    }
    if ($$ind{indid} eq '0') {
	$errstr = "$$dataset{pedigreefile}, $line $$dataset{pedlineno}: illegal person ID '0'";
	return undef;
    }
    if (($$ind{dadid} eq '0') != ($$ind{momid} eq '0')) {
	$errstr = "pedid $$ind{pedid}, person $$ind{indid} parents must both be either known or unknown";
	return (undef);
    }
    if ($$ind{indid} eq $$ind{dadid} || $$ind{indid} eq $$ind{momid}) {
        $errstr = "pedid $$ind{pedid}, person $$ind{indid} is coded as their own parent";
        return (undef);
    }
    if ($$ind{dadid} ne "0" && $$ind{dadid} eq $$ind{momid}) {
        $errstr = "pedid $$ind{pedid}, person $$ind{indid} parents are coded as the same person";
        return (undef);
    }
    
    foreach $trait (@{$$dataset{traitorder}}) {
	$traitcol = $$dataset{traits}{$trait}{col};
	push (@{$$ind{traits}}, $arr[$traitcol]);
	($$dataset{traits}{$trait}{flag} ne 'C' && defined ($$dataset{undefpheno})
	 && $arr[$traitcol] ne $$dataset{undefpheno})
	    and $$ind{phenotyped}++;
    }	      

    foreach $marker (@{$$dataset{markerorder}}) {
	($markercol, $idx) = @{$$dataset{markers}{$marker}}{qw/col idx/};
	if ($$ind{phased} && substr ($arr[$markercol], -1, 1, '') ne '|') {
	    $errstr = "$$dataset{pedigreefile}, line $$dataset{pedlineno}: individual $$ind{indid}, phased genotypes inconsistantly coded at marker $marker";
	    return (undef);
	}
	if (($arr[$markercol] eq '0') != ($arr[$markercol+1] eq '0')) {
	    $errstr = "$$dataset{pedigreefile}, line $$dataset{pedlineno}: individual $$ind{indid} is half-genotyped at marker $marker";
	    return (undef);
	}
	if ($arr[$markercol] ne '0') {
	    if (exists ($$dataset{markers}{$marker}{alleles})) {
		map {
		    if (!exists ($$dataset{markers}{$marker}{alleles}{$_})) {
			if ($$dataset{unkallelesok}) {
			    $$dataset{markers}{$marker}{alleles}{$_} = $KelvinDataset::ROUNDING_ERROR / 10;
			    warn ("WARNING, adding allele $_ for marker $marker\n");
			} else {
			    $errstr = "$$dataset{pedigreefile}, line $$dataset{pedlineno}: individual $$ind{indid} has unknown allele $_ for $marker";
			    return (undef);
			}
		    }
		} @arr[$markercol, $markercol+1];
	    }
	    $$ind{genotyped}++;
            $arr[$markercol] ne $arr[$markercol+1] and $$ind{heterozygous}++;
	}
	# Assign genotypes to the correct index in the individual's array of genotypes
	$$ind{markers}[$idx] = [ $arr[$markercol], $arr[$markercol+1] ];
    }	      
    return (bless ($ind, $class));
}

sub new_from_count
{
    my ($class, $dataset, $pedid, $traits, $labelmap, $aref) = @_;
    my $ind = { pedid => $pedid, indid => 1, dadid => 0, momid => 0,
		firstchildid => 0, patsibid => 0, matsibid => 0,
		origpedid => $pedid, origindid => 1, sex => undef,
                proband => 0, traits => [], markers => [], genotyped => 0,
                phenotyped => 0, heterozygous => 0, phased => 0,
                dataset => $dataset, makeped => 1 };
    my $trait;
    my $marker;
    my ($allele1, $allele2);
    my $va;

    $$ind{sex} = $$aref[2];
    foreach $trait (@{$$dataset{traitorder}}) {
	push (@{$$ind{traits}}, 'x');
	for ($va = 0; $va < scalar (@$traits); $va++) {
	    ($trait eq $$traits[$va]) and $$ind{traits}[-1] = $$aref[1][$va];
	}
    }
    ($allele1, $allele2) = split (//, $$aref[0]);
    for ($va = 0; $va < scalar (@{$$dataset{markerorder}}); $va++) {
	if (exists ($$labelmap[$va]{$allele1}) && exists ($$labelmap[$va]{$allele2})) {
	    $$ind{markers}[$va] = [@{$$labelmap[$va]}{($allele1, $allele2)}];
	    $$ind{genotyped}++;
            $$ind{markers}[$va][0] ne $$ind{markers}[$va][1] and $$ind{heterozygous}++;
	} else {
	    $$ind{markers}[$va] = [0, 0];
	}
    }
    return (bless ($ind, $class));
}

sub new_dummy
{
    my ($class, $dataset, $pedid, $indid, $sex) = @_;
    my $ind = { pedid => $pedid, indid => $indid, dadid => 0, momid => 0,
		firstchildid => 0, patsibid => 0, matsibid => 0,
		origpedid => $pedid, origindid => $indid, sex => $sex,
                proband => 0, traits => [], markers => [], genotyped => 0,
                phenotyped => 0, heterozygous => 0, phased => 0,
                dataset => $dataset, makeped => 'pre' };
    my $trait;
    my $marker;
    my ($allele1, $allele2);
    my $va;

    foreach $trait (@{$$dataset{traitorder}}) {
	push (@{$$ind{traits}}, 'x');
    }
    for ($va = 0; $va < scalar (@{$$dataset{markerorder}}); $va++) {
        $$ind{markers}[$va] = [0, 0];
    }
    return (bless ($ind, $class));
}

sub map
{
    my ($self, $newset) = @_;
    my $oldset = $$self{dataset};
    my $trait;
    my $marker;
    my $new = {traits => [], markers => [], genotyped => 0, phenotyped => 0,
               heterozygous => 0, dataset => $newset};
    
    map {
	$$new{$_} = $$self{$_};
    } qw/pedid indid dadid momid firstchildid patsibid matsibid origpedid
	origindid sex proband phased makeped/;

    foreach $trait (@{$$newset{traitorder}}) {
	if (exists ($$oldset{traits}{$trait})) {
	    push (@{$$new{traits}},  $$self{traits}[$$oldset{traits}{$trait}{idx}]);
	    ($$newset{traits}{$trait}{flag} ne 'C' && $$new{traits}[-1] != 0)
		and $$new{phenotyped}++;
	} else {
	    push (@{$$new{traits}}, 'x');
	}
    }
    foreach $marker (@{$$newset{markerorder}}) {
	if (exists ($$oldset{markers}{$marker})) {
	    push (@{$$new{markers}}, [ @{$$self{markers}[$$oldset{markers}{$marker}{idx}]} ]);
	    if ($$new{markers}[-1][0] ne '0') {
                $$new{genotyped}++;
                $$new{markers}[-1][0] ne $$new{markers}[-1][1] and $$new{heterozygous}++;
            }
	} else {
	    push (@{$$new{markers}}, [ 'x', 'x' ]);
	}
    }
    return (bless ($new, ref ($self)));
}

sub write
{
    my ($self) = @_;
    my $dataset = $$self{dataset};
    my $genosep = ($$self{phased} ? " | " : " ");

    if (! (defined ($$dataset{pedfh}) || $$dataset{writing})) {
	$dataset->writePedigreefile
	    or return (undef);
    }

    if ($$dataset{pedwriteformat} eq 'post') {
	# Write post-MAKEPED format. Don't care what the input format was.

	$$dataset{pedfh}->print (join (' ', @$self{qw/pedid indid dadid momid firstchildid patsibid matsibid sex proband/}, @{$$self{traits}}, map { "$$_[0]$genosep$$_[1]" } @{$$self{markers}}), "  Ped: $$self{origpedid}  Per: $$self{origindid}\n");

    } else {
	# Write pre-MAKEPED format from post-MAKEPED input.

	$$dataset{pedfh}->print (join (' ', @$self{qw/pedid indid dadid momid sex/}, @{$$self{traits}}, map { "$$_[0]$genosep$$_[1]" } @{$$self{markers}}), "\n");

    }
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
    if ($$dataset{traits}{$trait}{flag} ne 'C') {
	($$self{traits}[$idx] =~ /^(?:0|x)$/) and $$self{phenotyped}++;
	($value =~ /^(?:0|x)$/) and $$self{phenotyped}--;
    }
    $$self{traits}[$idx] = $value;
    return (1);
}

# Traits are scalar values, so returns the value directly
sub getTrait
{
    my ($self, $trait) = @_;
    my $dataset = $$self{dataset};
    my $idx;

    unless (exists ($$dataset{traits}{$trait})) {
	$errstr = "no trait $trait in dataset";
	return (undef);
    }
    $idx = $$dataset{traits}{$trait}{idx};
    return ($$self{traits}[$idx]);
}


sub setGenotype
{
    my ($self, $marker, $aref) = @_;
    my $dataset = $$self{dataset};
    my $idx;

    unless (exists ($$dataset{markers}{$marker})) {
	$errstr = "no marker $marker in dataset";
	return (undef)
    }
    $idx = $$dataset{markers}{$marker}{idx};
    @{$$self{markers}[$idx]} = @$aref;
    return (1);
}

# Genotypes are stored in array refs, so copy the values into a new array ref and
# return that
sub getGenotype
{
    my ($self, $marker) = @_;
    my $dataset = $$self{dataset};
    my $idx;
    my $aref = [];
    
    unless (exists ($$dataset{markers}{$marker})) {
	$errstr = "no marker $marker in dataset";
	return (undef)
    }
    $idx = $$dataset{markers}{$marker}{idx};
    @$aref = @{$$self{markers}[$idx]};
    return ($aref);
}

sub getAllGenotypes {
    # Returns *all* our markers.
    my ($self) = @_;
    
    return ($$self{markers});
}
sub setAllGenotypes {
    # Given an arrayref with marker info, sets all our markers en masse.
    # Used primarily when assembling MC-MC fully informative pedigree samples.
    my ($self, $markers) = @_;
    my $dataset = $$self{dataset};
    
    # minor sanity check - verify that we have the same number of markers as is
    # in the dataset
    my $datasetcount = scalar(@{$dataset->markerOrder()});
    unless (scalar(@$markers) == $datasetcount) {
        $errstr = "count mismatch between provided and dataset markers";
        return (undef)
    }
    
    $$self{markers} = $markers;
    return (1);
}

sub dataset
{
    my ($self) = @_;

    return ($$self{dataset});
}

sub pedid
{
    my ($self) = @_;

    return ($$self{pedid});
}

sub setPedid
{
    my ($self, $pedid) = @_;

    $$self{pedid} = $pedid;
    return (1);
}

sub indid
{
    my ($self) = @_;

    return ($$self{indid});
}

sub setIndid
{
    my ($self, $indid) = @_;

    $$self{indid} = $indid;
    return (1);
}

sub origindid
{
    my ($self) = @_;

    return ($$self{origindid});
}

sub dadid
{
    my ($self) = @_;

    return ($$self{dadid});
}

sub setDadid
{
    my ($self, $dadid) = @_;

    $$self{dadid} = $dadid;
    return (1);
}

sub momid
{
    my ($self) = @_;

    return ($$self{momid});
}

sub setMomid
{
    my ($self, $momid) = @_;

    $$self{momid} = $momid;
    return (1);
}

sub sex
{
    my ($self) = @_;

    return ($$self{sex});
}

sub setSex
{
    my ($self, $sex) = @_;

    $$self{sex} = $sex;
    return (1);
}

sub phenotyped
{
    my ($self) = @_;

    return ($$self{phenotyped});
}

sub genotyped
{
    my ($self) = @_;

    return ($$self{genotyped});
}

sub heterozygous
{
    my ($self) = @_;

    return ($$self{heterozygous});
}

sub phased
{
    my ($self) = @_;

    return ($$self{phased});
}

sub makeped
{
    my ($self) = @_;

    return ($$self{makeped});
}

sub founder
{
    my ($self) = @_;

    return ($$self{dadid} eq '0' && $$self{momid} eq '0');
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

1;
