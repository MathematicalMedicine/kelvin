#!/bin/bash

# manginl 04-27-07

# Shell script to run kelvin (either multipoint or two-point)
# in /tmp to avoid unnecessary network traffic. This script
# checks for appropriate files to be present, copies them to
# a tmp dir in /tmp, runs kelvin and copies the output back
# to the calling dir

# Checking that commandline args present (name of temp dir
# to create in /tmp
if [ -n "$1" ]; then :
else
	echo "This script requires a tmpdir name....exiting"
	exit
fi

# Creating temp dir
tmp_dir=$1
mkdir /tmp/$tmp_dir
final_dir=`pwd`

# Checking that all input files are present, if not, err msg
# then undo any action taken and exit
for i in kelvin.conf datafile.dat mapfile.dat markers.dat pedpost.dat; do
	if [ -e "$i" ]; then
		cp $i /tmp/$tmp_dir
	else
		echo "$i file not found...exiting"
		rm -rf /tmp/$tmp_dir
		exit
	fi
done

# Changing current working dir to /tmp/$tmp_dir
cd /tmp/$tmp_dir

# Starting kelvin run
date >kelvin.start
kelvin-0.0.29 kelvin.conf >kelvin.out </dev/null >&kelvin.err
date >kelvin.stop

# Cleaning up after kelvin run finishes

# copy output back to calling dir
cp * $final_dir/.

# change back to calling dir
cd $final_dir

# remove tmp_dir and contents
rm -rf /tmp/$tmp_dir

exit
