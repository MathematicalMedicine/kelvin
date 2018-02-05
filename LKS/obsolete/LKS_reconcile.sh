#!/bin/bash -e

# Reconcile LKS active Servers with their processes throughout the cluster. Update those that
# are missing so we can CleanOrphans.  This code is specific to the BCMM cluster as it uses
# the beo_exec utility, and beo_exec seems to require node names that cannot be intuited from
# local host name.
#
# Grab STUDY directive for database parameters. Use client.conf in current directory unless
# a different configuration file is specified as a parameter.
#
echo ' $Id$'

if test ! -z "$1"; then
    study=$(grep -i ^Study $1)
else
    study=$(grep -i ^Study client.conf)
fi
set -- $study

echo "set +e" > reconcile.sh
while read line
do
  echo -n "$i "
  i=$((i+1))
  set -- $line
  left=${1##Levi-Montalcini}
  host=${left%%.ccri.net}
  host=${host%%.crii.org}
  echo -n test -z \"\`beo_exec --node=n >>reconcile.sh || true
  echo -n $host ps h $2\`\" >>reconcile.sh
  echo -n ' && echo Update Servers set ExitStatus = 42 where ServerId = ' >>reconcile.sh
  echo $3 \\\; >>reconcile.sh
done <<<"$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select HostName, ProcessId, ServerId from Servers where StudyId = $2 AND ExitStatus IS NULL;")"
echo Running generated script...
source ./reconcile.sh >reconcile.sql
if test -s reconcile.sql; then
    echo Setting ExitStatus for `cat reconcile.sql | wc -l` servers
    #result=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="Source reconcile.sql;")"
else
    echo No missing servers detected
fi
#rm reconcile.sh reconcile.sql
