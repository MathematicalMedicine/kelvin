#!/bin/bash -eu

die() {
    echo $*
    exit -1
}

set -x

alias qrsh="qrsh -now no -cwd "

# Setup database tables
perl $KELVIN_ROOT/InitStudy.pl client.conf || die "InitStudy.pl client.conf failed"
perl $KELVIN_ROOT/InitStudy.pl server.conf || die "InitStudy.pl server.conf failed"

# Initial full run of client
qrsh $QUEUE -b y $KELVIN_ROOT/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0 || die "qrsh $KELVIN_ROOT/kelvin-2.2.0-study client.conf failed"

# Grab STUDY directive for database parameters
study=$(grep -i ^Study client.conf)
set -- $study

# Initial set of "new" trait positions is the original set so we can detect if _no_ splits ever occurred
cp client.conf client-newTP.conf
while :
do

  qsub $QUEUE -cwd -j y -t 1-9:1 -l virtual_free=2G -sync y -b y perl -I$KELVIN_ROOT $KELVIN_ROOT/merlin_lk_server.pl server.conf || die "qsub $KELVIN_ROOT/merlin_lk_server.pl server.conf failed"
  

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
  qrsh $QUEUE -b y $KELVIN_ROOT/kelvin-2.2.0-study client-newTP.conf --ProgressLevel 2 --ProgressDelaySeconds 0 || die "qsub $KELVIN_ROOT/kelvin-2.2.0-study client-newTP.conf failed"
  grep WARNING br.out || { break; }
  # Get the new set of trait positions
  TPs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select distinct RefTraitPosCM from Regions where StudyId = $2 AND RefTraitPosCM >= 0 AND PendingLikelihoods > 0;" | tr "\n" " ")
  grep -vi TraitPosition client.conf > client-newTP.conf
  # Add them one per line to be safe with line length
  for tp in $TPs;  do   echo "In the loop";  echo "TraitPosition $tp" >> client-newTP.conf; done
done

# Don't bother with a second run if there were no splits at all
diff -q client.conf client-newTP.conf && exit 0

qrsh $QUEUE -b y $KELVIN_ROOT/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0 || die "qrsh $KELVIN_ROOT/kelvin-2.2.0-study client.conf failed"

