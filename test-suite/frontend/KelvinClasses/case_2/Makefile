test :
	- rm -rf *.new*
	perl -I $(KELVIN_ROOT) classtest.pl
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.post.new.1 ped.post-test1 0 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.post.new.2 ped.post-test2 0 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.pre.new.1 ped.pre-test1 0 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.pre.new.2 ped.pre-test2 0 || log "Target $@ FAILED"
	- rm -rf *.new*

clean :
	- rm -rf *.new*
