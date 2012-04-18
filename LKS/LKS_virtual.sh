#!/bin/bash -eu
# Start with qrsh -l jobs=8,virtual_free=32G
if test -d /dev/shm/MySQL ; then
  echo /dev/shm/MySQL already exists!
  exit 1
fi
mkdir /dev/shm/MySQL
mkdir /dev/shm/MySQL/etc
mkdir /dev/shm/MySQL/data
mkdir /dev/shm/MySQL/data/mysqllog
cd ~/kit/mysql-5.1.50
make install
cp /usr/local/src/bcmmtools/LKS_virtual.cnf /dev/shm/MySQL/etc/my.cnf
export PATH=/dev/shm/MySQL/bin:/dev/shm/MySQL/libexec:$PATH
mysql_install_db --basedir=/dev/shm/MySQL --datadir=/dev/shm/MySQL/data
mysqld_safe --defaults-file=/dev/shm/MySQL/etc/my.cnf &
sleep 15
mysql --host=127.0.0.1 --user=root --execute="source /usr/local/src/bcmmtools/LKS_virtual.sql"
mysql --host=127.0.0.1 --user=adminLKS --pass=barfewww --database=LKS --execute="source /usr/local/src/bcmmtools/LKS_setup_tables.sql"
mysql --host=127.0.0.1 --user=adminLKS --pass=barfewww --database=LKS --execute="source /usr/local/src/bcmmtools/LKS_setup_trigger_proc.sql"
mysql --host=127.0.0.1 --user=adminLKS --pass=barfewww --database=LKS --execute="source /usr/local/src/bcmmtools/ModelParts.sql"
# Allow an hour for all tests to complete
#sleep 3600
sleep 300

