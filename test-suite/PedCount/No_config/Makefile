test :
	- rm -rf PC1_* kelvin.log
	$(KELVIN_ROOT)/PedCount.pl acid_test.post -post -count -write
	perl $(KELVIN_ROOT)/utils/wordDiff.pl PC1_Counts.Dat acid_test.counts-baseline 0 || log "Target $@ FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl PC1_Pedigrees.Dat acid_test.post-baseline 0 || log "Target $@ FAILED"

test-USE_DL :

clean :
	- rm -rf PC1_* kelvin.log
