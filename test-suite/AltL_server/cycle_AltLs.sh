perl ~/kelvin/trunk/InitStudy.pl $1client.conf
perl ~/kelvin/trunk/InitStudy.pl $1server.conf
while :
do
  ~/kelvin/trunk/kelvin-2.2.0 $1client.conf
  grep WARNING br.out
  if test $? -ne 0 ; then
    break
  fi
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  nq "~/kelvin/trunk/kelvin-2.2.0 $1server.conf"
  ~/kelvin/trunk/kelvin-2.2.0 $1server.conf
done
