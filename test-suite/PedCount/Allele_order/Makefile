test :
	- rm -rf O1_* NO1_*
	$(KELVIN_ROOT)/PedCount.pl -config ordered.conf -count -write=O
	$(KELVIN_ROOT)/PedCount.pl -config not-ordered.conf -count -write=NO
	perl $(KELVIN_ROOT)/utils/wordDiff.pl O1_Pedigrees.Dat NO1_Pedigrees.Dat .000001 || log "wordDiff for target $@ in $(PWD) FAILED"
	perl $(KELVIN_ROOT)/utils/wordDiff.pl O1_Counts.Dat NO1_Counts.Dat .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

test-USE_DL :

clean :
	- rm -rf O1_* NO1_*
