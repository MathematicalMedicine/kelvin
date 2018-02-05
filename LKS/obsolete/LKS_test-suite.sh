#!/bin/bash
# Start the virtual server parsing-out the job number information (e.g. Your job 
# 1342877 ("LKS_test-suite.sh") has been submitted), then get the host
# for that job once it starts running. Modify all conf files to point to the proper
# server instance, and run 'em while the database is up and available.

# Start the database server
qsub_output=$(qsub -l virtual_free=32G /usr/local/src/bcmmtools/LKS_test-server.sh)
echo qsub_output is $qsub_output
job_number=$(expr match "$qsub_output" 'Your job \([0-9]*\)')
echo job_number is $job_number

# Wait until the job actually starts running or 15 minutes passes
for i in {1..15}; do
    sleep 60;
    job_status=$(qstat | grep $job_number | tr -s ' ' | cut -d' ' -f 5)
    echo job_status is $job_status
    if test "$job_status" != "r" ; then
	if test -z "$job_status" ; then
	    echo "Job $job_number is gone - check the log. Aborting."
	    exit 1
	fi
	echo -n $i
	continue;
    fi
    queue_name=$(qstat | grep $job_number | tr -s ' ' | cut -d' ' -f 8);
    echo queue_name is $queue_name
    server_host=$(expr match "$queue_name" 'general@Levi-Montalcini\([0-9]*\)')
    echo server_host is $server_host
    break
done
if test -z "$server_host" ; then
    echo "It has been $i minutes and the job does not seem to be running. Aborting."
    qdel $job_number
    exit 1
fi
echo Job is running on queue $queue_name on n$server_host

# Get the standard MySQL parameters
study=$(grep -i ^Study /usr/local/src/bcmmtools/template-client.conf)
set -- $study

# Now give the job up to 5 minutes to get the instance running
for i in {1..5}; do
    sleep 60;
    Studies=$(mysql --host n$server_host --user $6 --password=$7 $5 --batch --skip-column-names --execute="Select count(*) from Studies;")
    echo Studies is $Studies
    if test -n "$Studies" ; then
	break;
    fi
    echo -n $i
done
if test "$Studies" -ne "0" ; then
    echo "Got the wrong value for Studies ($Studies)"
fi

# Now run all the tests, preferably in parallel
cd ~/kelvin/trunk/test-suite/likelihood_server
nq "cd SA_DT; make test"
nq "cd LC3; make test"

