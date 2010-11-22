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
set -x

# Don't just quit if nothing is available -- wait for it.
alias qrsh="qrsh -now no "
shopt -s expand_aliases

# These are for nodes other than Levi-Montalcini, where SGE is not available
if test "$HOSTNAME" != "Levi-Montalcini" ; then
    alias qrsh="bash -c "
    alias nq="echo Not submitting: "
fi

# Do the initialization only if there was no command line parameter
if test -z "$1" ; then
    # Setup database tables
    perl ~/kelvin/trunk/InitStudy.pl client.conf
    perl ~/kelvin/trunk/InitStudy.pl server-dataset1.conf
    perl ~/kelvin/trunk/InitStudy.pl server-dataset11.01.conf

    # Initial full run of client
    qrsh "cd `pwd`; ~/kelvin/trunk/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"

    # Initial set of "new" trait positions is the original set so we can detect if _no_ splits ever occurred
    cp client.conf client-newTP.conf
fi

# Grab STUDY directive for database parameters
study=$(grep -i ^Study client.conf)
set -- $study

# Setup the Single-Model RunTimes so bucket loading can be intelligent
SMRTs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Update PedigreePositions a, SingleModelRuntimes b set a.SingleModelEstimate = b.SingleModelRuntime, a.SingleModelRuntime = b.SingleModelRuntime where a.StudyId = $2 AND a.StudyId = b.StudyId AND a.PedigreeSId = b.PedigreeSId AND a.PedTraitPosCM = b.PedTraitPosCM;")
if test $SMRTs -ne 0 ; then
    # Assuming SOME finished, any that we missed should be treated as if they took too long..
    Singles=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Update PedigreePositions set SingleModelRuntime = 999999 where StudyId = $2 AND PedTraitPosCM <> 9999.99 AND SingleModelRuntime IS NULL;")
fi

while :
do
  # Enqueue a few servers and...
  nq "~/bcmmtools/run_server.sh server-dataset1"
  nq "~/bcmmtools/run_server.sh server-dataset1"
  nq "~/bcmmtools/run_server.sh server-dataset1"
  nq "~/bcmmtools/run_server.sh server-dataset1"
  nq "~/bcmmtools/run_server.sh server-dataset11.01"
  nq "~/bcmmtools/run_server.sh server-dataset11.01"
  nq "~/bcmmtools/run_server.sh server-dataset11.01"

  # Run single ones blocking further processing until most work is done
  qrsh "cd `pwd`; ~/bcmmtools/run_server.sh server-dataset1"
  qrsh "cd `pwd`; ~/bcmmtools/run_server.sh server-dataset11.01"
  # Make sure that nothing remains undone
  while :
  do
    ToDos=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select sum(PendingLikelihoods) from Regions where StudyId = $2;")
    if test $ToDos -eq 0 ; then
        break;
    fi
    Servers=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Servers where StudyId = $2 AND ExitStatus IS NULL;")
    if test $Servers -eq 0 ; then
	echo There are still Regions with incomplete work and no more servers are running!
	exit 1
    fi
    echo Waiting for servers to finish
    sleep 300
  done
  # Run the client to see if any splits occur
  qrsh "cd `pwd`; ~/kelvin/trunk/kelvin-2.2.0-study client-newTP.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
  grep WARNING br.out || { break; }
  # Get the new set of trait positions
  TPs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select distinct RefTraitPosCM from Regions where StudyId = $2 AND RefTraitPosCM >= 0 AND PendingLikelihoods > 0;" | tr "\n" " ")
  grep -vi TraitPosition client.conf > client-newTP.conf
  # Add them one per line to be safe with line length
  for tp in $TPs;  do   echo "In the loop";  echo "TraitPosition $tp" >> client-newTP.conf; done
done
# Don't bother with a second run if there were no splits at all
diff client.conf client-newTP.conf || qrsh "cd `pwd`; ~/kelvin/trunk/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"

# Get baselines of separate runs
nq "~/kelvin/trunk/kelvin-2.2.0-normal dataset1.conf --MapFile reference.map --BayesRatioFile br.reference-dataset1"
nq "~/kelvin/trunk/kelvin-2.2.0-normal dataset11.01.conf --MapFile reference.map --BayesRatioFile br.reference-dataset11.01"
