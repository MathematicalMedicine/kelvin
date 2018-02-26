#!/bin/bash -e
#
perl $KELVIN_ROOT/LKS/InitStudy.pl client.conf
perl $KELVIN_ROOT/LKS/InitStudy.pl server.conf

# Initial full run of client
qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"

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

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Run a single partial server - will do a subset of the pedigrees for all positions
qrsh "cd `pwd`; $KELVIN_ROOT/LKS/run_server.sh server-just-a-few"
ToDos=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Regions a, RegionModels b where a.AnalysisId = $AnalysisId AND a.RegionId = b.RegionId;")
echo Original client still has $ToDos models pending.

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Run the full client
rm study_*.dat || true
qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
grep WARNING br.out || true
echo Original client cannot complete analysis due to pending models.

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Try the reduced client
qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client-only-1.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
grep WARNING br.out || true
echo Reduced client finishes fine.

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Verify that positions are pending for the full client
TPs=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select distinct a.RefTraitPosCM from Regions a, RegionModels b where a.AnalysisId = $AnalysisId AND a.RegionId = b.RegionId AND a.RefTraitPosCM > 0.0;" | tr "\n" " ")
echo Positions $TPs are pending.

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Finish the pending full client models
qrsh "cd `pwd`; $KELVIN_ROOT/LKS/run_server.sh server"
ToDos=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Regions a, RegionModels b where a.AnalysisId = $AnalysisId AND a.RegionId = b.RegionId;")
echo Original client still has $ToDos models pending.

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Run the full client again to add more models (we know it splits)
rm study_*.dat || true
qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
grep WARNING br.out || true
echo Original client cannot complete analysis due to pending models.

# Reveal work to be done
FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")

# Run the reduced client to show that it is unaffected by the newly-added models
rm study_*.dat || true
qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client-only-1.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
grep WARNING br.out || true
echo Reduced client finishes fine.

# Run the full server and client a few more times to get final results
while :
do
  # Reveal work to be done
  FixFree=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call Q('FixFree');")
  qrsh "cd `pwd`; $KELVIN_ROOT/LKS/run_server.sh server"
  qrsh "cd `pwd`; $KELVIN_ROOT/kelvin-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0"
  grep WARNING br.out || { break; }
done

