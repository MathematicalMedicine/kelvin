#!/bin/sh
set -e

if [ -z "$CONF" ] ; then
    echo "no config file specified"
    exit -1
fi
if [ ! "$SGE_TASK_ID" -a "$SGE_TASK_LAST" ] ; then
    echo "no SGE variables in environment"
    exit -1
fi

MOD=`perl -we 'use strict;
               use KelvinConfig;
	       my $c = KelvinConfig->new ($ENV{CONF})
  	           or die ("KelvinConfig new failed, $KelvinConfig::errstr \n");
	       my $a = $c->isConfigured ("MODFile") or exit (0);
	       print ("$$a[0]\n");'`

BR=`perl -we 'use strict;
               use KelvinConfig;
	       my $c = KelvinConfig->new ($ENV{CONF})
  	           or die ("KelvinConfig new failed, $KelvinConfig::errstr \n");
	       print (${$c->isConfigured ("BayesRatioFile")}[0],"\n");'`

TP=`perl -we 'use strict;
	      use POSIX qw(ceil);
              use KelvinConfig;
              use KelvinDataset;
              my ($c, $d, $ar, $as, $s, $e, $i, $p, $x, %p, @p);
              $c = KelvinConfig->new ($ENV{CONF})
                  or die ("KelvinConfig new failed, $KelvinConfig::errstr \n");
	      $ar = $c->isConfigured ("TraitPositions")
                  or die ("no TraitPositions directive in $ENV{CONF}\n");
	      foreach $as (@$ar) {
                  if (($s, $e, $i) =
                      ($as =~ /(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?|end):(\d+(?:\.\d+)?)/)) {
                      if ($e eq "end") {
                          $d = KelvinDataset->new ({mapfile => ${$c->isConfigured ("MapFile")}[0]})
	                      or die ("KelvinDataset new failed, $KelvinDataset::errstr\n");
                          $e = ceil (${$d->getMarker (${$d->mapOrder}[-1])}{avgpos});
		          while ($e % $i) { $e++; } 
                      }
                      for ( $p = $s, $x = 1 ; $p <= $e; $p = $s + ($i * $x++)) { $p{$p} = ""; }
                  } else {
                      map { $p{$_} = ""; } split (/[, ]+/, $as);
                  }
              }
	      @p = sort { $a <=> $b } (keys (%p));
	      $i = int (scalar (@p) / $ENV{SGE_TASK_LAST}) + 1;
	      $s = $i * ($ENV{SGE_TASK_ID} - 1);
              print (join (",", splice (@p, $s, $i)), "\n");'`

if [ -z "$TP" ] ; then
    # This is not necessarily an error: when the number of tasks approaches the 
    # the number of trait positions, the last task may end up with no positions, 
    # depending on the math.
    echo "no trait positions to calculate, exiting"
    if [ $SGE_TASK_ID -eq $SGE_TASK_LAST ] ; then
	# Leave an empty BR file (and MOD file, if called for), so the merge doesn't croak
        echo "" > $BR.$SGE_TASK_ID
	[ "$MOD" ] && echo "" > $MOD.$SGE_TASK_ID
        exit 0
    fi
    exit -1
fi

set -x

echo "TraitPositions $TP" > $CONF.$SGE_TASK_ID
echo "BayesRatioFile $BR.$SGE_TASK_ID" >> $CONF.$SGE_TASK_ID
[ "$MOD" ] && echo "MODFile $MOD.$SGE_TASK_ID" >> $CONF.$SGE_TASK_ID
egrep -iv 'TraitPositions|MODFile|BayesRatioFile' $CONF >> $CONF.$SGE_TASK_ID
	       
$KELVIN_STUDY $CONF.$SGE_TASK_ID --ProgressLevel 2 --ProgressDelaySeconds 0

rm -f $CONF.$SGE_TASK_ID
