#!/usr/bin/perl -w
use strict;
use POSIX ":sys_wait_h";

$| = 0;

my $pid;
my $line;
my $vmsize;
my $vmpeak = 0;
my @ring = ();

my $loopdelay = 15;
my $ringsize = 100;

(($pid = fork ()) == -1)
    and die ("fork failed, $!\n");

if ($pid == 0) {
    exec (@ARGV) or die ("exec failed, $!\n");
}

while (1) {
    (waitpid ($pid, WNOHANG) == -1) and last;
    open (FP, "/proc/$pid/status") or last;
    while ($line = <FP>) {
	(($vmsize) = ($line =~ /vmsize:\s*(\d+)/i))
	    and last;
    }
    close (FP);
    ($vmpeak < $vmsize) and $vmpeak = $vmsize;
    (scalar (@ring) >= $ringsize)
	and shift (@ring);
    push (@ring, $vmsize);
    sleep ($loopdelay);
}

print (map { "$_\n" } ("Last $ringsize values:", @ring));
print ("Peak: $vmpeak\n");

exit (0);
