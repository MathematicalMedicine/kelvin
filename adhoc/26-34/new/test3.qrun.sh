#!/bin/bash
#$ -N test3
#$ -S /bin/bash
#$ -o test3.stdo.out
#$ -e test3.stderr.out
#$ -cwd
CDIR=`pwd`
. $HOME/.bash_profile
TMPD=/tmp/$LOGNAME/$JOB_ID
this_node=`/bin/uname -n`
echo "TMPD is $this_node:$TMPD"
mkdir -p $TMPD

cp test3.ped $TMPD
cp test3.mrk $TMPD
cp test3.dat $TMPD
cp test3.run $TMPD
cp test3.map $TMPD
cd $TMPD

/home/lbrz/bin/kelvin test3.run

cp $TMPD/test3.lods.out $CDIR
cp $TMPD/test3.avghet.out $CDIR
cp $TMPD/test3.ppl.out $CDIR
cd /tmp
rm -rf $TMPD
