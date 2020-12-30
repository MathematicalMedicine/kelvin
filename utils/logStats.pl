#!/usr/bin/perl
#
# Extracts performance statistics from kelvin logs in the normal test-suite hierarchy. Produces
# SQL statements to insert them in MySQL.
#
# Matches any of the following date/time formats from the "time" command (different platforms!):
#
# 1h2m3.4s
# 1m2.3s
# 1.2s
# 1:2:3.4s
# 1:2.3s
#
# Expects to be invoked as:
#
# find . \( -name kelvin\*.log -or -name cycle.out \) -exec egrep -i -e "[kK]elvin.* V[0-9\.]+ edit [0-9\.]+|^real\s+[0-9]+"  \{\} /dev/null \; | \
#   perl -S -n logStats.pl >>Inserts.sql
#
#print "For $_:";
if (/(kelvin.*)\.log:(\d{2}\/\d{2}\/\d{2} \d{2}\:\d{2}\:\d{2}) kelvin.* V(.*) edit (.*) built .*$/) {
#    print "K1=yes\n";
    ($RunDate, $Kelvin, $Build) = ($2, $1.'-'.$3, $4);
} else {
#    print "K1=no\n";
}
if (/cycle\..*out:(\d{2}\/\d{2}\/\d{2} \d{2}\:\d{2}\:\d{2}) kelvin.* V(.*) edit (.*) built .*$/) {
#    print "C1=yes\n";
    ($RunDate, $Kelvin, $Build) = ($1, 'kelvin-study-'.$2, $3);
} else {
#    print "K1=no\n";
}
if (/(dynamic|fixed)-grid\/(.*)\/kelvin.*log:real\s+(((\d*)[h:])?(\d*)[m:])?(\d+\.\d+)(s)?$/) {
#    print "K2=yes\n";
    print "Insert into testRuns (Kelvin, Build, Host, Test, Realtime, RunDate) values ('$Kelvin', '$Build', '".$ENV{HOSTNAME}."','$1/$2', ".(($5*3600)+($6*60)+$7).", STR_TO_DATE('$RunDate', '%y/%m/%d %H:%i:%s'));\n";
} else {
#    print "K2=no\n";
}
if (/(likelihood_server)\/(.*)\/cycle.*out:real\s+(((\d*)[h:])?(\d*)[m:])?(\d+\.\d+)(s)?$/) {
#   print "C2=yes\n";
    print "Insert into testRuns (Kelvin, Build, Host, Test, Realtime, RunDate) values ('$Kelvin', '$Build', '".$ENV{HOSTNAME}."','$1/$2', ".(($5*3600)+($6*60)+$7).", STR_TO_DATE('$RunDate', '%y/%m/%d %H:%i:%s'));\n";
} else {
#    print "C2=no\n";
}
    
