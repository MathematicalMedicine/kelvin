#!/usr/bin/env perl
use strict;
use warnings;
use IO::File;

#
# IO::File::Kelvin: A convenience wrapper around IO::File that only returns
# nonblank, noncomment lines, and counts all lines read (including blank and
# comment lines).
#
# Copyright (C) 2010, 2022 Mathematical Medicine LLC
# 
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
package IO::File::Kelvin;

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
        # Drop leading or trailing spaces, new lines, carriage returns and inline comments
        $line =~ s/(^\s*|(?:\s*\#.*)?\s*[\n\r]+$)//g;
	if ($line) {
	    # This is clumsy. How to do it better?
	    (ref ($ref) eq 'SCALAR' && $$ref =~ /^\d+$/) and $$ref += $lineno;
	    return ($line);
	}
    }
    return (undef);
}
