#!/bin/bash -eu

set -x

# These are for nodes other than Levi-Montalcini
if test "$HOSTNAME" != "Levi-Montalcini" ; then
    shopt -s expand_aliases
    alias qrsh="bash -c"
    alias nq="echo Not submitting: "
fi

# Do the normal 2pt run for comparison in a normal queue
nq "~/kelvin/trunk/kelvin-2.2.0-normal client-normal-2pt.conf --ProgressLevel 2 --ProgressDelaySeconds 0"

# Setup database tables
perl ~/kelvin/trunk/InitStudy.pl client.conf
perl ~/kelvin/trunk/InitStudy.pl server-mp.conf
perl ~/kelvin/trunk/InitStudy.pl server-2pt-init.conf

# Initial full run of client
qrsh "cd `pwd`; ~/kelvin/trunk/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"

# Grab STUDY directive for database parameters
study=$(grep -i ^Study client.conf)
set -- $study

# Set all 2pt PedigreePositions to one marker
echo "Update PedigreePositions set MarkerCount = 1 where StudyId = $2 AND PedigreeSId in (180);" | mysql --host=$4 --user=$6 --password=$7 $5

# Initial set of "new" trait positions is the original set
cp client.conf client-newTP.conf
while :
do
  # Enqueue a few servers and...
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-mp"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  nq "~/bcmmtools/run_server.sh server-2pt-run"
  # Run single blocking ones to prevent further processing until all work is done
  qrsh "cd `pwd`; ~/bcmmtools/run_server.sh server-mp"
  qrsh "cd `pwd`; ~/bcmmtools/run_server.sh server-2pt-run"
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
  grep WARNING br.out-all-pedigrees || { break; }
  # Get the new set of trait positions
  TPs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select distinct RefTraitPosCM from Regions where StudyId = $2 AND RefTraitPosCM >= 0 AND PendingLikelihoods > 0;" | tr "\n" " ")
  grep -vi TraitPosition client.conf > client-newTP.conf
  # Add them one per line to be safe with line length
  for tp in $TPs;  do   echo "In the loop";  echo "TraitPosition $tp" >> client-newTP.conf; done
done
# Don't bother with a second run if there were no splits at all
diff client.conf client-newTP.conf || qrsh "cd `pwd`; ~/kelvin/trunk/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
# Run a 2pt for just the 2pt pedigrees for comparison purposes
nq "~/kelvin/trunk/kelvin-2.2.0-study client-just-2pt.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
