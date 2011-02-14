#!/bin/bash -eu
/usr/local/src/kelvin-2.3.1/kelvin-2.3.1-study $1.conf --ProgressLevel 2 --ProgressDelaySeconds 0 &
study=$(grep -i ^Study $1.conf)
set -- $study
kid=$!
wait $kid
kidStatus=$?
mysql --host $4 --user $6 --password=$7 $5 --execute="Update Servers set StopTime = Now(), ExitStatus = $kidStatus where HostName = '`hostname`' AND ProcessId = $kid;"