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
# $Id$
#
shopt -s expand_aliases
set -x

# Any queuing modifier, like using -q johntest
qmods=""

# Handle unique job submission situations.
case $HOSTNAME in
    Levi-Montalcini )
    # Don't just quit if nothing is available -- wait for it.
    alias qrsh="qrsh -now no $qmods"
    ;;
    opt* ) # OSC's Opteron cluster
    ;;
    * ) # Everything else
    alias qrsh="bash -c "
    alias nq="bash -c "
    ;;
esac

# Do the initialization only if there was no command line parameter
if test -z "$1" ; then
    # Setup database tables
    perl $KELVIN_ROOT/InitStudy.pl client.conf
    perl $KELVIN_ROOT/InitStudy.pl server.conf

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

if test -z "$StudyId"; then
    echo "ERROR - STUDY directive in configuration file specified StudyLabel ($2) that was not found in the database, exiting!"
    exit 2
fi

# Setup the Single-Model RunTimes so bucket loading can be intelligent
SMRTs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Update PedigreePositions a, SingleModelRuntimes b set a.SingleModelEstimate = b.SingleModelRuntime, a.SingleModelRuntime = b.SingleModelRuntime where a.StudyId = $StudyId AND a.StudyId = b.StudyId AND a.PedigreeSId = b.PedigreeSId AND a.PedTraitPosCM = b.PedTraitPosCM;")
if test "$SMRTs" != "0" ; then
    # Assuming SOME finished, any that we missed should be treated as if they took too long..
    Singles=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Update PedigreePositions set SingleModelRuntime = 999999 where StudyId = $StudyId AND PedTraitPosCM <> 9999.99 AND SingleModelRuntime IS NULL;")
fi

while :
do
  # Enqueue no more servers than DB server threads until we're sure they're needed (and then by hand)
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &
  nq "$KELVIN_ROOT/run_server.sh server $qmods" &

  # Run a single one blocking further processing until most work is done
  qrsh "cd `pwd`; $KELVIN_ROOT/run_server.sh server"
  # Make sure that nothing remains undone
  while :
  do
    ToDos=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select sum(PendingLikelihoods) from Regions where StudyId = $StudyId;")
    if test $ToDos -eq 0 ; then
        break;
    fi
    Servers=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Servers where StudyId = $StudyId AND ExitStatus IS NULL;")
    if test $Servers -eq 0 ; then
	echo There are still Regions with incomplete work and no more servers are running!
	exit 1
    fi
    echo Waiting for servers to finish
    sleep 300
  done
  # Run the client to see if any splits occur
  qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client-newTP.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
  grep WARNING br.out || { break; }
  # Get the new set of trait positions
  TPs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select distinct RefTraitPosCM from Regions where StudyId = $StudyId AND RefTraitPosCM >= 0 AND PendingLikelihoods > 0;" | tr "\n" " ")
  grep -vi TraitPosition client.conf > client-newTP.conf
  # Add them one per line to be safe with line length
  for tp in $TPs;  do   echo "In the loop";  echo "TraitPosition $tp" >> client-newTP.conf; done
done
# Don't bother with a second run if there were no splits at all
diff client.conf client-newTP.conf || qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
