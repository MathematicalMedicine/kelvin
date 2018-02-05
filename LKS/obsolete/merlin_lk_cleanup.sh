#!/bin/sh

# The typical name for cycle script output
rm -f cycle.out
# The client configuration with trait positions, created between client and server runs
rm -f client-newTP.conf
# Output files created by the merlin_lk_server array job
rm -f perl.o*
# Output from the kelvin client runs
rm -f kelvin-*-study.o*
# Output from split_client_job and split_client_merge runs
rm -f split_client_job.sh.o* split_client_merge.sh.o*
# Bayes ratio and conf files from busted split_client_job runs
rm -f br.out.[1234] mod.out.[1234] client.conf.[1234]
# Intermediate files from client runs
rm -f study_*.dat
