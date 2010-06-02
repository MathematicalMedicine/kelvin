#!/usr/bin/perl -w 
use strict;
use IO::File;

#
# IO::File::Kelvin: A convenience wrapper around IO::File that only returns
# nonblank, noncomment lines, and counts all lines read (including blank and
# comment lines).
#
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
	$line =~ s/^\s*([^\#\n]*)\#?.*$/$1/s;
	if ($line) {
	    # This is clumsy. How to do it better?
	    (ref ($ref) eq 'SCALAR' && $$ref =~ /^\d+$/) and $$ref += $lineno;
	    return ($line);
	}
    }
    return (undef);
}
