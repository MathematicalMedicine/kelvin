test :
	- rm -f *.full *.cc
	(time $(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post-full --BayesRatioFile br.full --PPLFile ppl.full --MODFile mod.full ) >kelvin.full 2>&1
	(time $(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post-cc --BayesRatioFile br.cc --PPLFile ppl.cc --MODFile mod.cc ) >kelvin.cc 2>&1
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.full br.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.cc br.full .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.full ppl.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.cc ppl.full .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl mod.full mod.out-baseline .000001
	perl $(KELVIN_ROOT)/utils/wordDiff.pl mod.cc mod.full .000001

test-USE_DL :
	echo Compilation not appropriate for this test.

test-FIXED :
	echo Fixed-grid comparison not available

clean :
	- rm -f *.full *.cc
