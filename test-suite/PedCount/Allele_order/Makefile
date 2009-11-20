test :
	- rm -rf O1_* NO1_*
	$(KELVIN_ROOT)/PedCount.pl -config ordered.conf -count -write=O
	$(KELVIN_ROOT)/PedCount.pl -config not-ordered.conf -count -write=NO
	diff O1_Pedigrees.Dat NO1_Pedigrees.Dat
	diff O1_Counts.Dat NO1_Counts.Dat

test-USE_DL :

clean :
	- rm -rf O1_* NO1_*
