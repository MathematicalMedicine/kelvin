#!/bin/sh
set -e

if [ -z "$CONF" ] ; then
    echo "no configuration file specified"
    exit -1
fi

BR=`perl -we 'use strict;
               use KelvinConfig;
               my ($c, $a);
               $c = KelvinConfig->new ($ENV{CONF})
                   or die ("KelvinConfig new failed, $KelvinConfig::errstr \n");
               $a = $c->isConfigured ("BayesRatioFile") or exit (0);
               print ("$$a[0]\n");'`

MOD=`perl -we 'use strict;
               use KelvinConfig;
               my ($c, $a);
	       $c = KelvinConfig->new ($ENV{CONF})
  	           or die ("KelvinConfig new failed, $KelvinConfig::errstr \n");
	       $a = $c->isConfigured ("MODFile") or exit (0);
	       print ("$$a[0]\n");'`

cp $BR.1 $BR
[ "$MOD" ] && cp $MOD.1 $MOD
TASK_ID=2
while [ $TASK_ID -le $COUNT ] ; do
    egrep -v 'Version|Chr' $BR.$TASK_ID >> $BR
    [ "$MOD" ] && egrep -v 'Version|Chr' $MOD.$TASK_ID >> $MOD
    TASK_ID=`expr $TASK_ID + 1`
done

TASK_ID=1
while [ $TASK_ID -le $COUNT ] ; do
    rm -f $BR.$TASK_ID $MOD.$TASK_ID
    TASK_ID=`expr $TASK_ID + 1`
done
