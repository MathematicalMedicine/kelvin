for i in $(seq 1 1000); do qsub -cwd -l virtual_free=1G -b n -N TQLC1P${i} TPCT_sT_poly.sh; done
