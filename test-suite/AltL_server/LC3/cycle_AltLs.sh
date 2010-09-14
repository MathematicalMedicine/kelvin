perl ~/kelvin/trunk/InitStudy.pl client.conf
perl ~/kelvin/trunk/InitStudy.pl server.conf
while :
do
  ~/kelvin/trunk/kelvin-2.2.0-study client.conf
  grep WARNING br.out
  if test $? -ne 0 ; then
    break
  fi
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0-study server.conf"
  ~/kelvin/trunk/kelvin-2.2.0-study server.conf
done
