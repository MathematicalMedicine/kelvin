#!/bin/bash -eu
# Start with qrsh -l jobs=8,virtual_free=32G

MYSQL_BASE=/opt/mysql-5.5.12

if test -d /dev/shm/MySQL ; then
  echo /dev/shm/MySQL already exists!
  exit 1
fi
mkdir -p /dev/shm/MySQL/data
cp $MYSQL_BASE/support-files/ramdisk.cnf /dev/shm/MySQL/my.cnf
export PATH=$MYSQL_BASE/bin:$PATH

mysql_install_db --basedir=$MYSQL_BASE --datadir=/dev/shm/MySQL/data --user=$LOGNAME
mysqld_safe --defaults-file=/dev/shm/MySQL/my.cnf &
sleep 15
mysql --host=localhost --user=root --execute="source /usr/local/src/bcmmtools/LKS_virtual.sql"
mysql --host=127.0.0.1 --user=adminLKS --pass=barfewww --database=LKS --execute="source /usr/local/src/bcmmtools/LKS_setup_tables.sql"
mysql --host=127.0.0.1 --user=adminLKS --pass=barfewww --database=LKS --execute="source /usr/local/src/bcmmtools/LKS_setup_trigger_proc.sql"
mysql --host=localhost --user=adminLKS --pass=barfewww --database=LKS --execute="source /usr/local/src/bcmmtools/ModelParts.sql"
