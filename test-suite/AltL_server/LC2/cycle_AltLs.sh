perl ~/kelvin/trunk/InitStudy.pl client.conf
perl ~/kelvin/trunk/InitStudy.pl server.conf
while :
do
  ~/kelvin/trunk/kelvin-2.2.0-study client.conf
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
  ~/kelvin/trunk/kelvin-2.2.0-study server.conf
done
