#!/bin/bash -e
#
# Run the likelihood server version of kelvin as cleverly as possible to minimize the overall
# runtime. Study-specific versions are common, with changes reflecting the various servers
# that have to run to perform all the work between client passes.
#
# Really intended to run on the head node so that batch jobs can be submitted, but it will run
# (slowly) on other nodes for simple tests.
#
# Initially run without any parameters, as it gets everything it needs from the STUDY directive
# line in the client.conf file. If a parameter is specified (no matter what it is), then the
# initialization steps will be skipped under the assumption that this is a resumed run.
#
# It is definitely best to have run SMRT.sh first to generate and load Single Model RunTimes
# for all pedigree/positions. Otherwise servers could be gobbling-up 50 models, each of which 
# takes days to run, and running them in serial.
#
# Copyright (C) 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
# $Id$
#
shopt -s expand_aliases
set -x

# Any queuing modifier, like using -q johntest
qmods=""

# Handle unique job submission situations.
case $HOSTNAME in
    levi-Montalcini* )
        # Don't just quit if nothing is available -- wait for it.
        alias qrsh="qrsh -now no $qmods "
        lks_server_count=${lks_server_count-8}
    ;;
    opt* ) # OSC's Opteron cluster
        lks_server_count=${lks_server_count-8}
    ;;
    master | node* ) # StarCluster
        alias qrsh="qrsh -now no $qmods "
        alias nq="qsub -cwd "
        lks_server_count=${lks_server_count-2}
    ;;
    * ) # Everything else
        alias qrsh="bash -c "
        alias nq="bash -c "
        lks_server_count=${lks_server_count-2}
    ;;
esac

# Do the initialization only if there was no command line parameter
if test -z "$1" ; then
    # Setup database tables
    perl $KELVIN_ROOT/LKS/InitStudy.pl client.conf
    perl $KELVIN_ROOT/LKS/InitStudy.pl server.conf

    # Initial full run of client
    qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"

    # Initial set of "new" trait positions is the original set so we can detect if _no_ splits ever occurred
    cp client.conf client-newTP.conf
fi

# Grab STUDY directive for database parameters
study=$(grep -i ^Study client.conf)
set -- $study

# Get the actual StudyId from the StudyLabel. This MUST succeed if InitStudy.pl is worth its salt.
StudyId=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="Select StudyId from Studies where StudyLabel = '$2';")

# Get the AnalysisId from the database.
AnalysisId=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="Select AnalysisId from Analyses where StudyId = $StudyId AND PedigreeRegEx = '$8' AND PedigreeNotRegEx = '$9'")

if test -z "$StudyId"; then
    echo "ERROR - STUDY directive in configuration file specified StudyLabel ($2) that was not found in the database, exiting!"
    exit 2
fi

c=1
while :
do
  # Reveal work to be done
    FixFree=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

  # Enqueue no more servers than DB server threads until we're sure they're needed (and then by hand)
  for ((servs=1; servs<lks_server_count; servs++)); do
    qrsh "cd `pwd`; $KELVIN_ROOT/LKS/run_server.sh server " & 
  done
  # Run a single one blocking further processing until most work is done
  qrsh "cd `pwd`; $KELVIN_ROOT/LKS/run_server.sh server"
  while :
  do
    # Reveal work to be done
    FixFree=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

    ToDos=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Regions a, RegionModels b where a.AnalysisId = $AnalysisId AND a.RegionId = b.RegionId;")
    if test $ToDos -eq 0 ; then
        break;
    fi
    Servers=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Servers where StudyId = $StudyId AND ExitStatus IS NULL;")
    if test $Servers -eq 0 ; then
        echo There are still Regions with incomplete work and no more servers are running!
        exit 1
    fi
    echo Waiting for servers to finish
  done
  # Reveal work to be done
  FixFree=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")
  # Run the client to see if any splits occur
  qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client-newTP.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
  cp br.out br.out.$c
  c=$((c+1))
  grep WARNING br.out || { break; }
  # Get the new set of trait positions
  TPs=$(mysql --host=$4 --user=$6 --password=$7 $5 --batch --skip-column-names --execute="Select distinct a.RefTraitPosCM from Regions a, RegionModels b where a.AnalysisId = $AnalysisId AND a.RegionId = b.RegionId AND a.RefTraitPosCM > 0.0;" | tr "\n" " ")
  grep -vi TraitPosition client.conf > client-newTP.conf
  # Add them one per line to be safe with line length
  for tp in $TPs;  do   echo "In the loop";  echo "TraitPosition $tp" >> client-newTP.conf; done
done
# Don't bother with a second run if there were no splits at all
diff client.conf client-newTP.conf || qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
