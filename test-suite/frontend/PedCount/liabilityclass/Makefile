test :
	- rm -f ped.post.counted ped.post.counts *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.log 2>&1
	#diff -ib ped.post.counted ped.post.counted-baseline
	#diff -ib ped.post.counts ped.post.counts-baseline
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-baseline .000001 || log "Target $@ FAILED"

clean :
	- rm -f ped.post.counted ped.post.counts *.out

