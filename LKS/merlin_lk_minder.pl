#!/usr/bin/perl
use strict;
use warnings FATAL => qw(all);
use DBI;
use POSIX ":sys_wait_h";

my $max_children = 4;
my $loop_delay = 5;
my $tooldir = '/export/local/bcmmtools/lks';
my $cleanupscript = "$tooldir/merlin_lk_cleanup.sh";
my $cyclescript = "$tooldir/merlin_lk_cycle.sh";

my %child_pids;
my @workdirs;
my $ret;

( -x $cleanupscript && -x $cyclescript )
    or die ("cycle scripts are missing not executable in '$tooldir'\n");

while (scalar (@ARGV)) {
    (-d $ARGV[0])
	or die ("'$ARGV[0]' is not a directory\n");
    (-f "$ARGV[0]/client.conf" && -f "$ARGV[0]/server.conf")
	or die ("'$ARGV[0]' doesn't contain basic configuration files\n");
    check_study ("$ARGV[0]/client.conf");
    push (@workdirs, shift (@ARGV));
}

while (scalar (@workdirs) > 0 || scalar (keys (%child_pids)) > 0) {
    while (scalar (keys (%child_pids)) < $max_children && scalar (@workdirs)) {
	spawn_child (\%child_pids, $workdirs[0]) and shift (@workdirs);
    }
    foreach (keys (%child_pids)) {
	$ret = waitpid ($_, WNOHANG);
	if ($ret == $_) {
	    print ("$_ : finished $child_pids{$_}, status $?\n");
	    delete ($child_pids{$_});
	} elsif ($ret == -1) {
	    print ("WARNING: waitpid $_ failed, $!\n");
	}
    }
    sleep ($loop_delay);
}


sub spawn_child
{
    my ($href, $dir) = @_;
    my $pid;
    my $ret;

    $pid = fork ();
    if (! defined ($pid)) {
	print ("WARNING: fork failed, $!\n");
	sleep (10);
    } elsif ($pid != 0) {
	print ("$pid : starting $dir\n");
	$$href{$pid} = $dir;
	return (1);
    }
    chdir ($dir)
	or die ("$$: chdir '$dir' failed, $!\n");

    $ret = system ("$cyclescript > cycle.out 2>&1");
    if ($ret == -1) {
	print ("$$: WARNING system first '$cyclescript' failed, $!\n");
	exit ($ret);
    } elsif ($ret != 0) {
	print ("$$: WARNING first '$cyclescript' exited in error\n");
	exit ($ret);
    }

    $ret = system ("$cleanupscript > /dev/null 2>&1");
    if ($ret == -1) {
	print ("$$: WARNING system '$cleanupscript' failed, $!\n");
	exit ($ret);
    } elsif ($ret != 0) {
	print ("$$: WARNING second '$cleanupscript' exited in error\n");
	exit ($ret);
    }

    $ret = system ("$cyclescript > cycle.out 2>&1");
    if ($ret == -1) {
	print ("$$: WARNING system second '$cyclescript' failed, $!\n");
	exit ($ret);
    } elsif ($ret != 0) {
	print ("$$: WARNING second '$cyclescript' exited in error\n");
	exit ($ret);
    }


    exit (0);
}


sub check_study
{
    my ($conf) = @_;
    my $line;
    my @arr;
    my $dbh;
    my $aref;

    open (FILE, $conf) or die ("open '$conf' failed, $!\n");
    while ($line = <FILE>) {
	$line =~ /^Study/ and last;
    }
    ($line) or die ("'$conf' doesn't contain a Study directive\n");
    @arr = split (' ', $line);
    (scalar (@arr) == 9 && $arr[0] =~ /study/i)
	or die ("badly formed Study directive in '$conf'\n");
    ($arr[2] =~ /client/i)
	or die ("'$conf' is a client config, but says '$arr[2]' in the Study directive\n");
    $dbh = DBI->connect ("DBI:mysql:hostname=$arr[3];database=$arr[4]", $arr[5], $arr[6])
	or die ("DBI connect to $arr[3], db $arr[4] as $arr[5] failed, $DBI::errstr\n");

#   InitStudy runs from inside Merlin cycle script, so the study won't be in the DB yet
#    $aref = $dbh->selectall_arrayref ("select StudyId from Studies where StudyLabel = '$arr[1]'")
#	or die ("DBI selectall for study '$arr[1]' failed, $DBI::errstr\n");
#    (scalar (@$aref) != 1)
#	and die ("Study '$arr[1]' does not exist in database '$arr[4]'\n");

    return (1);
}
