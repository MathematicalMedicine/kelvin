perl ~/kelvin/trunk/InitStudy.pl nonpoly-client.conf
perl ~/kelvin/trunk/InitStudy.pl nonpoly-server.conf
while :
do
  ~/kelvin/trunk/kelvin-2.2.0-study nonpoly-client.conf
  grep WARNING br.out
  if test $? -ne 0 ; then
    break
  fi
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  nq "~/bcmmtools/run_server.sh nonpoly-server"
  ~/kelvin/trunk/kelvin-2.2.0-study nonpoly-server.conf
done
