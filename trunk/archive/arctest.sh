export OMP_NUM_THREADS=1
rm MAKE-FAILED* | true
make -j 4 -f ../Makefile.common
for i in MAKE-RUNNING-*; do mv $i ${i/RUNNING/FAILED}; done
