#!/usr/bin/env perl

# Adds an allocatable "database" resource to our compute nodes. Necessary for
# Kelvin-LKS to function properly!

package LKS_setupSGEDB v0.1.0;

use v5.20;
use warnings 'FATAL';
use lib ($ENV{TOOLPATH} ? split (/:/, $ENV{TOOLPATH}) : ' ');
use English;

use Getopt::Long qw(:config posix_default auto_version no_ignore_case);
use List::MoreUtils qw(any);
use File::Temp;

use Explorer qw(import)
our @EXPORT_OK = qw(init_lks_db find_in_path);

$OUTPUT_AUTOFLUSH = 1;
my $svnid = '$Id$';


# Global - need to know if any qconf/qhost/etc. commands should be run using
# 'sudo'
my $sudo_enabled = 0;


sub _add_master_resource {
    # Adds the "database" complex attribute to our grid (if not already
    # present) so that it can be referenced by compute nodes and allocated by
    # jobs later.
    
    print("adding master 'database' attribute to SGE...\n");
    my $sgeconf = [];
    open(my $qconf_sc, "-|", _sudo("qconf -sc"))
            || die("Could not run qconf: $!");
    while (<$qconf_sc>) {
        chomp;
        if ($ARG =~ /^database/) {
            print("...'database' attribute already exists; step skipped.\n");
            return
        } else {
            push(@$sgeconf, $ARG);
        }
    }
    close($qconf_sc);
    push(@$sgeconf, "database db INT <= YES YES 0 0");
    push(@$sgeconf, "", "");
            # qconf -Mc insists on an extra empty newline at the end of its
            # input and will error out if that line is not present (sigh)
    
    _load_tempfile(_sudo("qconf -Mc"), join("\n", @$sgeconf));
    print("...'database' attribute added.\n");
}


sub _add_node_database_resource {
    # Modifies the compute node named so that it references our "database"
    # complex attribute, so that the database facilities on said node can be
    # allocated and tracked.
    my ($nodename) = @ARG;
    
    print("adding node $nodename as a database node...\n");
    my $nodeconf = {};
    open(my $qconf_se, "-|", _sudo("qconf -se $nodename"))
            || die("Could not run qconf: $!");
    
    # we need to convert our output into what amounts to a hash of hashes so
    # that we can reliably search it for the parameter we want and append it if
    # needed.
    # the output is a tad messy - `qconf -se $node > foo; qconf -Me foo`,
    # despite this ostensibly being The Way To go, would actually fail! it
    # chokes on "continuation lines" (multiline strings ending in backslashes)
    # in the output, and some of the key/value pairs that are given. real nice
    # Kwality Kontrol we've got here.
    my $linetemp = "";
    while (<$qconf_se>) {
        chomp;
        if ($ARG =~ / \\$/) {
            # continuation line (or the last part of a continuation line)
            # just tack it on for now; we'll remove characters later
            $linetemp = $linetemp . "_NEWLINE_" . $ARG;
        } elsif ($linetemp) {
            # end of a multiline string
            $linetemp = $linetemp . "_NEWLINE_" . $ARG;
            # get rid of excess spaces and the backslashes, then add.
            $linetemp =~ s/ *\\*_NEWLINE_ *//;
            _hashify_nodeconf($nodeconf, $linetemp);
            $linetemp = "";
        } else {
            _hashify_nodeconf($nodeconf, $ARG);
        }
    }
    close($qconf_se);
    
    unless (any { $ARG eq "database=1" } @{$$nodeconf{"complex_values"}}) {
        # add the database parameter we want
        if ($$nodeconf{"complex_values"}[0] eq "NONE") {
            pop(@{$$nodeconf{"complex_values"}});
        }
        push(@{$$nodeconf{"complex_values"}}, "database=1");
        
        # reassemble our hash into a string
        my $revisedconf = "";
        while (my ($key, $value) = each(%$nodeconf)) {
            $revisedconf = $revisedconf . "$key  " . join(",", @$value) . "\n";
        }
        
        _load_tempfile(_sudo("qconf -Me"), $revisedconf);
        print("...$nodename now marked as a database node.\n");
    } else {
        print("...$nodename already marked as a database node; skipping.\n");
    }
}

sub _hashify_nodeconf {
    # Given a hashref we've been adding key/value pairs to, and a string to be
    # added to the hashref, convert that string into a key/value pair and
    # append it to the hashref.
    my ($nodeconf, $newval) = @ARG;
    
    # check the new value and make sure it's not one of the parameters that
    # qconf -Me chokes on (no, really, qconf seriously will hand you parameters
    # that it doesn't take back, because reasons)
    # we have to do this here because some of those verboten parameters are
    # multiline strings, so it has to take place after those have been joined
    # back together. ist das nicht wunderfubar?
    return if ($newval =~ /^load_values/);
    return if ($newval =~ /^processors/);
    
    # we're good, so let's divide and append it
    my ($key, $value) = split(/\s+/, $newval);
    $$nodeconf{$key} = [split(/,/, $value)];
    
    return;
}

sub _load_tempfile {
    # Given two strings - one a command that expects a filename, and one
    # representing the contents of a file to be fed to that command, creates a
    # temporary file and feeds it to that command.
    my ($command, $contents) = @ARG;
    
    my ($loadfile, $loadfilename) = File::Temp::tempfile(CLEANUP => 1);
    print $loadfile $contents;
    system("$command $loadfilename");
    ($? != 0) && die(sprintf("$command failed with signal %d", $? >> 8));
    close($loadfile);
}

sub _get_host_info {
    # Get our list of nodes and information regarding same from qhost, and
    # return it as a hash.
    
    my $hostinfo = {};
    open(my $qhost, "-|", _sudo("qhost")) || die("could not run qhost: $!");
    while (<$qhost>) {
        chomp;
        next if /^HOSTNAME/;
        next if /^--------/;
        next if /^global  /;
        my $hostline = [split];
        my $hostname = shift(@$hostline);
        $$hostinfo{$hostname} = $hostline;
    }
    close($qhost);
    return ($hostinfo);
}

sub _sudo {
    # Given a string that represents a command, if the --use_sudo option is
    # turned on, returns that string prefixed to use sudo.
    my ($command) = @ARG;
    
    if ($sudo_enabled) {
        return "sudo " . $command;
    } else {
        return $command;
    }
}

sub setup_db_resource {
    # Given an arrayref of nodes to set up as database nodes, adds them to the
    # scheduler as same (and also configures the scheduler, if needed).
    my ($nodelist, $skip_add_master, $use_sudo) = @ARG;
    
    if ($use_sudo) { $sudo_enabled = 1; }
    
    print("Checking for possible runtime problems...\n");
    open(my $testrun, "-|", _sudo("qconf -help"))
            || die("qconf test run failed: $!");
    
    my $hostinfo = _get_host_info();
    foreach my $node (@$nodelist) {
        unless ($$hostinfo{$node}) { die("node $node not found in host list"); }
    }
    print("...we should be good to go!\n");
    
    # no problems that we know of, so let's go ahead and make this happen.
    unless ($skip_add_master) { _add_master_resource(); }
    foreach my $node (@$nodelist) { _add_node_database_resource($node); }
    print("All done!\nRemember to make sure your MySQL distribution is present on each DB node, in the same directory name on each one!\n\n");
}

sub _print_help_and_exit {
    # Prints help and exits with the given status code.
    my ($exitstatus) = @ARG;
    
    print <<END_OF_HELP;

usage: LKS_setupSGEDB.pm [--host <hostname>] [--host <hostname>] [...] [--use_sudo] [--skip_add_master]


Adds an allocatable "database" resource to the designated compute nodes.
Necessary for Kelvin-LKS to function properly!

The script does NOT actually copy the required MySQL distribution to each node;
that has to be done by you, by hand.

This script should be run from a node that can run the `qhost` and `qconf`
commands. In most cases, this will be your head node.


    --host <hostname> - Add the "database" resource to the specified SGE
                        node. More than one may be specified on the command
                        line.

    --use_sudo - Runs `qhost` and `qconf` with sudo instead of directly. Useful
                 if you don't particularly feel like running this script as a
                 privileged user.

    --skip_add_master - Skip attempting to create the "database" resource in
                        SGE. Saves some time if you've already done this
                        previously (will not cause damage if you omit it for
                        already-configured queue masters though).

END_OF_HELP
    
    exit($exitstatus);
}

sub main {
    my ($skip_add_master, $use_sudo, $help);
    my $hostlist = [];
    Getopt::Long::GetOptions(
        "host=s" => $hostlist,
        "skip_add_master" => \$skip_add_master,
        "use_sudo" => \$use_sudo,
        "help" => \$help,
    ) || _print_help_and_exit(1);
    
    ($help) && _print_help_and_exit(0);
    unless ($$hostlist[0]) {
        print("Need to specify at least one host!\n\n");
        _print_help_and_exit(1);
    }
    
    setup_db_resource($hostlist, $skip_add_master, $use_sudo);
}


main() unless caller(0);
1;
