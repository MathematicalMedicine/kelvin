#!/bin/bash
# Fire-up the virtual database server and leave it running by executing
# this script using: qsub -l virtual_free=32G /usr/local/src/bcmmtools/LKS_test-server.sh

# Safety check - don't want to wipe-out someone else's analysis...
if test -d /dev/shm/MySQL ; then
  echo /dev/shm/MySQL already exists!
  exit 1
fi
# Startup the database server
/usr/local/src/bcmmtools/LKS_virtual.sh
# Give time for it to be used, say about 8 hours
#sleep 28800
sleep 1800
# Shutdown the database server
mysqladmin --host=127.0.0.1 --user=adminLKS --pass=barfewww shutdown
sleep 10
# Clean-out the directory
rm -rf /dev/shm/MySQL/
