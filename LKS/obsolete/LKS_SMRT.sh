#!/bin/bash -e

#
# This script performs a multipoint "Single-Model RunTimes" (SMRT)
# sizing run on a single chromosome for all pedigrees using the
# multipoint marker count provided on the command line. If you need 2pt
# or a mix, ask Bill for help setting up the run, as it is complicated.
#
# This script should be invoked in the directory with the configuration
# files, typically a chromosome subdirectory (e.g. chr6). It requires 
# a multipoint dynamic grid kelvin configuration file as the first 
# parameter (typically kelvin.conf), e.g.:
#
# /usr/local/bin/LKS_SMRT.sh ../kelvin.conf
#
# If the script is executed in a chrN subdirectory, it will use a common
# "hundred block" StudyId number for all sibling directories in order to
# make combined results analysis feasible.
#
# This script MUST BE EXECUTED ON THE HEADNODE as it
# submits batch jobs for each pedigree.
#
# Once all submitted jobs have completed, the script will collect
# runtime data from the log files and load it into MySQL on Walker in
# the SingleModelRuntimes table under StudyId.
#
# After the runtime data is loaded, stored procedure SMRT_stats is
# executed and the results stored in the current directory as
# SMRT_stats.out.  These are the statistics for the single run.
# You may want to invoke the procedure for all StudyIds in your
# analysis in order to get summary statistics.
#
# To execute copies of this script in all subdirectories of an
# analysis, go to the parent directory of the chrN subdirectories
# and issue a command like:
#
# for i in chr*; do cd $i; /usr/local/src/bcmmtools/LKS_SMRT.sh ../kelvin.conf >/dev/null & cd ..; done
#
# Note that standard output is redirected to /dev/null. This is to keep the
# noise level low, as otherwise the output of all subshells would be displayed
# at once. Output from the 'for' command itself and any errors will still be
# displayed. Once the run is complete, stay in the parent directory and issue
# the command:
#
# /usr/local/src/bcmmtools/LKS_SMRT_summary.sh kelvin.conf
#
# NOTICE that for this example, the reference to kelvin.conf is for the current 
# directory since the script is executing in the parent direction. This will
# give you a summary of statistics for all chrN subdirectories.
#
# If you really know what you are doing, you can specify a second
# parameter to use as a StudyId or "hundred block". Generally don't do this.
#

echo ' $Id: LKS_SMRT.sh 298 2012-03-28 17:12:25Z whv001 $'

if test "$HOSTNAME" != "Levi-Montalcini" ; then
    echo "FATAL - This must run on the headnode in order to submit the SMRT jobs." 1>&2
    exit 1
fi
if test -z "$1" ; then
    echo "FATAL - You must specify the multipoint dynamic grid kelvin configuration file (with relative path, if needed)." 1>&2
    exit 1
else
    ConfigFile=$1
fi
if ! test -e "$1" ; then
    echo "FATAL - Cannot find specified  multipoint dynamic grid kelvin configuration file ($1)." 1>&2
    exit 1
fi

# Allow for an optional second parameter to use as the StudyId or "hundred block"
OverrideNumber=$2

# Grab the standard template STUDY directive for database parameters
Study=$(grep -i ^Study /usr/local/src/bcmmtools/template-client.conf)
set -- $Study
Host=$4
Database=$5
Username=$6
Password=$7

FullDirectory=$(pwd)
BaseDirectory=$(basename $FullDirectory)
ParentDirectory=${FullDirectory%$BaseDirectory}
ChromosomeNumber=${BaseDirectory#chr}

if [[ "$ChromosomeNumber" =~ ^[0-9]+$ ]] ; then
    if test -z "$OverrideNumber"; then
        # Maybe create and then retrieve the "Hundred Block"
	CreateHundredBlock=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Insert into HundredBlockSizingRuns (Directory) value ('$ParentDirectory');" 2>/dev/null) || true
	HundredBlock=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Select HundredBlock*100 from HundredBlockSizingRuns where Directory = '$ParentDirectory';") || true
    else
	HundredBlock=$OverrideNumber
    fi

    if [[ "$HundredBlock" =~ ^[0-9]+00$ ]] ; then
	StudyId=$(($HundredBlock+$ChromosomeNumber))
	echo -n "Results for this collection of runs will be under StudyIds $HundredBlock"
        echo -n "01 thru $HundredBlock"
	echo "23."
    else
	echo "FATAL ($BaseDirectory instance) - Failed to acquire sane Hundred Block." 1>&2
	exit 1
    fi
else
    if test -z "$OverrideNumber"; then
        # Maybe create and then retrieve the StudyId
	CreateSingleStudyId=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Insert into SingleSizingRuns (Directory) value ('$FullDirectory');") || true
	StudyId=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Select StudyId from SingleSizingRuns where Directory = '$FullDirectory';") || true
    else
	StudyId=$OverrideNumber
	echo "The results for this run will be under StudyId $StudyId."
    fi

    if ! [[ "$StudyId" =~ ^[0-9]+$ ]] ; then
	echo "FATAL ($BaseDirectory instance) - Failed to acquire sane StudyId [$StudyId]." 1>&2
	exit 1
    fi
fi

# Grab the number of multipoint markers
Multipoint=$(grep -i ^Multipoint $ConfigFile)
set -- $Multipoint
Markers=$2

if test -e SMRT_$Markers; then
    echo "Removing old SMRT_$Markers subdirectory."
    rm -rf SMRT_$Markers
fi

# See if this is an imprinting run
Imprinting=$(grep -i ^Imprinting $ConfigFile) || true
if test -n "$Imprinting"; then
  ImprintingFlag='y'
  ImprintingName='Imp'
else
  ImprintingFlag='n'
  ImprintingName='NoImp'
fi

# Get the liability class count
LiabilityClass=$(grep -i ^Liability $ConfigFile) || true
set -- $LiabilityClass
LiabilityClassCnt=$2
if test -z "$LiabilityClassCnt"; then
  LiabilityClassCnt=1
fi

# See if this is a DT or QT run
QT=$(grep -i ^QT $ConfigFile) || true
if test -n "$QT"; then
    TraitName='QT'
else
    TraitName='DT'
fi

# Add the Studies row ignoring failure
AddStudy=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Insert into Studies (StudyId, ReferenceMapId, LiabilityClassCnt, ImprintingFlag) values ($StudyId, -1, $LiabilityClassCnt, '$ImprintingFlag');" 2>/dev/null) || true

# Grab pedigree filename
PedigreeFile=$(grep -i ^Pedigree $ConfigFile) || true
if test -z "$PedigreeFile" ; then
    PedigreeFileName=pedfile.dat
else
    set -- $PedigreeFile
    PedigreeFileName=$2
fi

mkdir SMRT_$Markers

cat $ConfigFile /usr/local/src/bcmmtools/SMRT.$TraitName.$ImprintingName.conf-part | grep -vi Multipoint >SMRT_$Markers.conf
echo "Multipoint $Markers" >>SMRT_$Markers.conf

name=SMRT_$RANDOM

JobsRemaining=0
for i in $(cut -f1 -d\  $PedigreeFileName | sort | uniq)
do
    JobsRemaining=$((JobsRemaining+1))
    grep ^$i\  $PedigreeFileName >SMRT_$Markers/pedfile_$i.dat
    (nq "(time /usr/local/bin/kelvin SMRT_$Markers.conf --PedigreeFile SMRT_$Markers/pedfile_$i.dat --BayesRatioFile SMRT_$Markers/br.out-$i --ProgressDelaySeconds 0 --ProgressLevel 1) >SMRT_$Markers/kelvin.out-$i 2>&1" -N $name )>/dev/null
done

echo "Waiting while $JobsRemaining batch jobs named $name complete (starting at `date`.)"
echo ""
echo "Remember that if the runs are taking an inordinately long time to complete,"
echo "like over an hour, you can find most-recently-updated log files with the"
echo "command 'ls -lsrt $FullDirectory/SMRT_$Markers/kelvin.out-*'."
echo "You can then type those files to see the percent progress."
echo ""

# Think about -hold_jid $name here...would work, but we'd have to write a temp file of commands that would include the password, or re-parse it...hmmm.

while [ $JobsRemaining -gt 0 ]
do
    sleep 60
    JobsRemaining=$(qstat | grep $name | wc -l)
    echo -ne "\r$JobsRemaining..."
done

echo Batch jobs complete, collecting results and loading into database under StudyId $StudyId.
cd SMRT_$Markers
# See if there were any kelvin errors before we hose down the work surface
errors=$(egrep -e FATAL -e ERROR -e WARNING kelvin.out-[0-9]*) || true
if ! test -z "$errors" ; then
    echo "FATAL ($BaseDirectory instance) - ERROR or WARNING messages found in kelvin log files, check the logs in the SMRT_$Markers subdirectory" 1>&2
fi
perl /usr/local/src/bcmmtools/SMRT_extract.pl $StudyId 0 $Markers
RemoveOld=$(mysql --host $Host --user $Username --password=$Password $Database --batch --skip-column-names --execute="Delete from SingleModelRuntimes where StudyId = $StudyId AND MarkerCount = $Markers;")
mysql --host $Host --user $Username --password=$Password $Database --execute='source SMRT_Inserts.sql'

# Generate statistics for this run only
echo Load complete, generating statistics to `pwd`/SMRT_stats.out
mysql --host $Host --user $Username --password=$Password $Database --execute="call SMRT_stats($StudyId, $StudyId)" >SMRT_stats.out
cd ..
rm SMRT_[0-9]*.[oev][0-9]*
rm SMRT_$Markers/br.out-[0-9]*
rm SMRT_$Markers/pedfile_[0-9]*.dat
rm kelvin_[0-9]*_memory.dat 2>/dev/null
echo Done.
