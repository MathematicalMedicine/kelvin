#!/bin/bash -e

# This script performs summarizing queries on a set of
# "Single-Model Runtimes" (SMRT) results in the MySQL database
# on Walker.
#
# It is invoked with two parameters, the low and high
# StudyIds to include in the results.
#

#set -x

echo " $Id$"

if test -z "$2" ; then
    echo "You must specify low and high StudyIds to include in the summary."
    exit 1
fi

low=$1
high=$2

# Grab the standard template STUDY directive for database parameters
Study=$(grep -i ^Study ~/bcmmtools/template-client.conf)
set -- $Study
Host=$4
Database=$5
Username=$6
Password=$7

mysql --host $Host --user $Username --password=$Password $Database --execute="call SMRT_stats($low, $high);"



