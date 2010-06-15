test :
	- rm -f *.full *.cc
	(time $(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post-full --BayesRatioFile br.full --PPLFile ppl.full --MODFile mod.full ) >kelvin.full 2>&1
	(time $(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post-cc --BayesRatioFile br.cc --PPLFile ppl.cc --MODFile mod.cc ) >kelvin.cc 2>&1
	diff -ib -I '^# Version .*' br.out-baseline br.full
	diff -i br.full br.cc
	diff -ib -I '^# Version .*' ppl.out-baseline ppl.full
	diff -i ppl.full ppl.cc
	diff -ib -I '^# Version .*' mod.out-baseline mod.full
	diff -i mod.full mod.cc

clean :
	- rm -f *.full *.cc
