perl ~/kelvin/trunk/InitStudy.pl client.conf
perl ~/kelvin/trunk/InitStudy.pl server.conf
while :
do
  ~/kelvin/trunk/kelvin-2.2.0 client.conf
  grep WARNING br.out
  if test $? -ne 0 ; then
    break
  fi
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 server.conf"
  ~/kelvin/trunk/kelvin-2.2.0 server.conf
done
