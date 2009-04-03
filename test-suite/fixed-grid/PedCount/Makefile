test :
	- rm -rf PC1_* kelvin.out
	$(KELVIN_ROOT)/PedCount.pl acid_test.post -post -count -write
	diff -i PC1_Counts.Dat acid_test.counts-baseline
	diff -i PC1_Pedigrees.Dat acid_test.post-baseline

clean :
	- rm -rf PC1_* kelvin.out
