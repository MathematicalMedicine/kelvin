#!/bin/bash -e

#
#

echo ' $Id$'

if test -z "$1" ; then
    echo "FATAL - You must specify the multipoint dynamic grid kelvin configuration file (with relative path, if needed)." 1>&2
    exit 1
else
    ConfigFile=$1
fi

# Grab the standard template STUDY directive for database parameters
Study=$(grep -i ^Study /usr/local/src/bcmmtools/template-client.conf)
set -- $Study
Host=$4
Database=$5
Username=$6
Password=$7

FullDirectory=$(pwd)

# Find the original "Hundred Block"...
HundredBlock=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Select HundredBlock*100 from HundredBlockSizingRuns where Directory = '$FullDirectory/';") || true
if [[ "$HundredBlock" =~ ^[0-9]+00$ ]] ; then
    LowStudyId=$(($HundredBlock+1))
    HighStudyId=$(($HundredBlock+23))
else
    # ...or find the original StudyId
    StudyId=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Select StudyId from SingleSizingRuns where Directory = '$FullDirectory';") || true
    LowStudyId=$StudyId
    HighStudyId=$StudyId
    if ! [[ "$StudyId" =~ ^[0-9]+$ ]] ; then
	echo "FATAL - Cannot find sizing information for this directory." 1>&2
	exit 1
    fi
fi

echo Generating statistics to SMRT_stats.out
mysql --host $Host --user $Username --password=$Password $Database --execute="call SMRT_stats($LowStudyId, $HighStudyId)" >SMRT_stats.out

echo Done.
