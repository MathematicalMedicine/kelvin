#!/bin/bash -eu
perl ~/kelvin/trunk/InitStudy.pl client.conf
perl ~/kelvin/trunk/InitStudy.pl server.conf
while :
do
  ~/kelvin/trunk/kelvin-2.2.0-study client.conf --ProgressLevel 2 --ProgressDelaySeconds 0
  grep WARNING br.out
  if test $? -ne 0 ; then
    break
  fi
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  nq "~/bcmmtools/run_server.sh server"
  ~/bcmmtools/run_server.sh server
done
