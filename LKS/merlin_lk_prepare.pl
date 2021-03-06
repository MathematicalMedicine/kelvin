#!/usr/bin/env perl
# Copyright (C) 2011, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
use strict;
use FindBin;
use lib split(/:+/, "!:$ENV{'KELVIN_ROOT'}:NO_KELVIN_ROOT:$ENV{'TOOLPATH'}:$FindBin::Bin");
use warnings;
use KelvinConfig;
use DBI;

# Used globally
my $dsn;
my $dbh;

my $configfile = shift (@ARGV);
my $config;
my %study;
my $liability;
my $sth;
my $sql;
my $lc_select = '';
my $aref;
my $rows;
my $modelPartTable = 'DModelParts';
my $retries = 0;

# build insert ... select statement based on config (liability classes, mostly)
# execute
# disconnect from database

print (ts(), "$0 starting on $ENV{HOSTNAME} in $ENV{PWD}, pid $$". (exists ($ENV{JOB_ID}) ? ", job ID $ENV{JOB_ID}" : ""). "\n");

($configfile && -f $configfile)
    or die ("usage: $0 <configfile>\n");
print (ts(), "Using configfile $configfile\n");
$config = KelvinConfig->new ($configfile)
    or die ("KelvinConfig->new failed: $KelvinConfig::errstr\n");
$aref = $config->isConfigured ("Study")
    or die ("no Study directive in config file '$configfile'\n'");
@study{qw/label role dbhost dbname dbuser dbpasswd pedregex pednotregex/} = split (' ', $$aref[0]);
($study{role} =~ /server/i)
    or die ("Study role in config file '$configfile' is not 'server'\n");
($config->isConfigured ("QT") || $config->isConfigured ("QTT")) 
    and $modelPartTable = 'QModelParts';
$aref = $config->isConfigured ("LiabilityClasses");
$liability = defined ($aref) ? $$aref[0] : 1;

$dsn = "DBI:mysql:host=$study{dbhost};database=$study{dbname}";
$dbh = DBI->connect ($dsn, $study{dbuser}, $study{dbpasswd},
#		     {AutoCommit => 0})
		     {AutoCommit => 0, PrintError => 0})
    or die ("DBI connect to '$dsn' as $study{dbuser} failed, $DBI::errstr\n");

db_get_studyid ($dbh, \%study);  

print (ts(), "StudyLabel $study{label}, Id $study{id}");
print (($liability > 1) ? ", $liability liability classes\n" : ", 1 liability class\n");

db_cleanup_old ($dbh, \%study);

$lc_select = join (", ", map { "LC${_}MPId" } (1 .. $liability));

$sql = "insert into LGModels (LGModelId, StudyId, PedigreeRegEx, " .
       "  PedigreeNotRegEx, " . $lc_select. ") ".
       "select distinct 0, p.StudyId, ?, ?, " . $lc_select . " ".
       "from PedigreePositions p, Models m, ". $modelPartTable. " d ".
       "where p.StudyId = ? ".
       "  and p.PedigreeSId regexp ? ".
       "  and p.PedigreeSId not regexp ? ".
       "  and m.Likelihood is NULL ".
       "  and m.ServerId is NULL ".
       "  and p.PedPosId = m.PedPosId ".
       "  and m.LC1MPId = d.MPId ".
       "  and d.DGF != -1";

$sql =~ s/  */ /g;
print ("SQL is '$sql'\n");

($sth = $dbh->prepare ($sql))
    or die ("prepare insert into LGModels failed, $DBI::errstr\n");

print (ts(), "Inserting rows into LGModels\n");
while (1) {
    if (! ($rows = $sth->execute (@study{qw/pedregex pednotregex id pedregex pednotregex/}))) {
	($DBI::errstr !~ /try restarting transaction/)
	    and die ("insert into LGModels failed, $DBI::errstr\n");
	print ("retrying insert into LGModels\n");
        $retries++;
    } else {
	last;
    }
}
$dbh->commit or die ("commit insert into LGModels failed, $DBI::errstr\n");
print (ts(), "Inserted $rows into LGModels with $retries retries\n");

exit (0);


sub db_get_studyid
{
    my ($dbh, $study) = @_;
    my $sth;
    my $aref;

    (($sth = $dbh->prepare ("select StudyId from Studies ".
                            "where StudyLabel = ? "))
     && $sth->execute ($$study{label}))
        or die ("select from Studies failed, $DBI::errstr\n");
    (defined ($aref = $sth->fetchrow_arrayref))
        or die ("select from Studies returned no data\n");
    $$study{id} = $$aref[0];
    return (1);
}


sub db_cleanup_old 
{
    my ($dbh, $study) = @_;

    (($sth = $dbh->prepare ("delete from LGModels ".
                            "where StudyId = ? ".
                            "  and PedigreeRegEx = ? ".
                            "  and PedigreeNotRegEx = ?"))
     && ($rows = $sth->execute (@$study{qw(id pedregex pednotregex)})))
        or die ("delete from LGModels failed, $DBI::errstr\n");
    $dbh->commit or die ("commit delete from LGModels failed, $DBI::errstr\n");
    print (ts(), "Deleted $rows from LGModels\n");
    return (1);
}


sub ts
{
    my @arr = localtime ();
    return (sprintf ("%02d/%02d %02d:%02d:%02d ", $arr[4]+1, @arr[3,2,1,0]));
}
