#!/bin/bash
# Copyright (C) 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

~/kelvin/trunk/kelvin-study $1.conf --ProgressLevel 2 --ProgressDelaySeconds 0 &
study=$(grep -i ^Study $1.conf)
set -- $study
kid=$!
wait $kid
kidStatus=$?
mysql --host $4 --user $6 --password=$7 $5 --execute="Update Servers set StopTime = Now(), ExitStatus = $kidStatus where HostName = '`hostname`' AND ProcessId = $kid;"
