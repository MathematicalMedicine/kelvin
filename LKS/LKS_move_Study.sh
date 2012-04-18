#!/bin/bash -e
#
# Move an LKS study between databases, possibly on completely different servers.
#
# Execute in the directory for the Study to dump. References TWO FILES:
# client.conf and target.conf. The first is the "from" client configuration file
# and the second is the "to" client configuration file. Only the database access
# information from these two files is used. The easy way to use this is just copy
# your client.conf to target.conf and edit the database access details. This is
# NOT FOR MOVING STUDIES WITHIN A DATABASE. Use the stored procedure 
# MoveStudy(fromStudyId, toStudyId) for that.
#
# Will find the first StudyId that is available in BOTH databases. Attempts to 
# preserve the current "from" database StudyId, but if it is in conflict with
# a study in the "to" database, will find a StudyId that is currently unused in
# both databases and use it instead.  When this happens, THE "FROM" STUDYID
# IS CHANGED IN THE "FROM" DATABASE AS WELL in order to make the extract work.
#
# Once the "from" and "to" StudyIds are in agreement, all related rows are extracted
# to text files and loaded in the "to" database.  The only exception is that ALL
# MODEL PARTS rows (from both DModelParts and QModelParts) are moved.
#
# $Id$
#
# William Valentine-Cooper, 2012-03
#

shopt -s expand_aliases
set -x

# Get the two sets of database connect information.
# Grab STUDY directive for database parameters
from=$(grep -i ^Study client.conf)
set -- $from
fl=$2; fh=$4; fd=$5; fu=$6; fp=$7
to=$(grep -i ^Study target.conf)
set -- $to
tl=$2; th=$4; td=$5; tu=$6; tp=$7

# Get the "from" StudyId from the "from" StudyLabel.
fi=$(mysql --host=$fh --user=$fu --password=$fp $fd --batch --skip-column-names --execute="Select StudyId from Studies where StudyLabel = '$fl';")
# The "copy" StudyId starts out as the "from" StudyId
ti=$fi
# Loop until we get an acceptable StudyId in both databases
while :
do
  tl=$(mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="Select StudyLabel from Studies where StudyId = $ti;")
  if test -z "$tl"; then
      echo "Found a good target StudyId of $ti"
      break
  fi
  # Find the next available source target StudyId
  while :
  do
    ti=$((ti+1))
    tl=$(mysql --host=$fh --user=$fu --password=$fp $fd --batch --skip-column-names --execute="Select StudyLabel from Studies where StudyId = $ti;")
    if test -z "$tl"; then
	echo "Found another candidate target StudyId of $ti"
	break
    fi
  done
done

# If we had to change StudyIds, update everything in the source database.
if test "$fi" != "$ti"; then
  rows=$(mysql --host=$fh --user=$fu --password=$fp $fd --batch --skip-column-names --execute="call MoveStudy($fi,$ti);")
fi

# Extract metadata creations only so we can allow them to fail.  We'll keep processing going with --force
mysqldump  --host $fh --user $fu --password=$fp $fd --no-data --add-drop-table=false >LKS_metadata_to_StudyId_$ti.sql
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_metadata_to_StudyId_$ti.sql;" --force && true

# Now extract all data related to the target StudyId
mysqldump  --host $fh --user $fu --password=$fp $fd Studies --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Studies_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Analyses --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Analyses_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Maps --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Maps_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Markers --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Markers_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd MapMarkers --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_MapMarkers_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Pedigrees --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Pedigrees_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Positions --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Positions_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd PedigreePositions --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_PedigreePositions_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Regions --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Regions_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd Servers --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_Servers_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd LGModels --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_LGModels_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd SingleSizingRuns --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_SingleSizingRuns_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd SingleModelRuntimes --no-create-info --skip-lock-tables --skip-triggers --where="StudyId = $ti" >LKS_SingleModelRuntimes_to_StudyId_$ti.sql
# And everything indirectly related to the target StudyId
mysqldump  --host $fh --user $fu --password=$fp $fd Models --no-create-info --skip-lock-tables --skip-triggers --where="PedPosId in (Select PedPosId from PedigreePositions where StudyId = $ti)" >LKS_Models_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd MarkerSetLikelihood --no-create-info --skip-lock-tables --skip-triggers --where="PedPosId in (Select PedPosId from PedigreePositions where StudyId = $ti)" >LKS_MarkerSetLikelihood_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd RegionModels --no-create-info --skip-lock-tables --skip-triggers --where="RegionId in (Select RegionId from Regions where StudyId = $ti)" >LKS_RegionModels_to_StudyId_$ti.sql
mysqldump  --host $fh --user $fu --password=$fp $fd ServerPedigrees --no-create-info --skip-lock-tables --skip-triggers --where="ServerId in (Select ServerId from Servers where StudyId = $ti)" >LKS_ServerPedigrees_to_StudyId_$ti.sql
# And all model parts on general principles
mysql --host=$fh --user=$fu --password=$fp $fd --batch --skip-column-names --execute="call UnloadModelParts();" >LKS_ModelParts_to_StudyId_$ti.sql

exit

# Finally load everthing into the new database in the order that will not violate referential integrity.
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Studies_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Analyses_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Maps_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Markers_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_MapMarkers_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Pedigrees_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Positions_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_PedigreePositions_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Regions_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Servers_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_LGModels_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_SingleSizingRuns_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_SingleModelRuntimes_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_ModelParts_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_Models_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_MarkerSetLikelihood_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_RegionModels_to_StudyId_$ti.sql;"
mysql --host=$th --user=$tu --password=$tp $td --batch --skip-column-names --execute="source LKS_ServerPedigrees_to_StudyId_$ti.sql;"
