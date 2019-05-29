test :
	- rm -f br.out br.tmp br.split ppl.out kelvin.out*
	$(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	$(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post.1 --LocusFile locus.dat.1 \
		--BayesRatioFile br.split >>kelvin.out 2>&1
	$(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post.2 --LocusFile locus.dat.2 \
		--BayesRatioFile br.tmp >>kelvin.out 2>&1
	fgrep -v Version br.tmp | perl -pe 's/Seq: 1/Seq: 2/' >> br.split
	$(TEST_KELVIN) kelvin.conf --PedigreeFile ped.post.3 --LocusFile locus.dat.3 \
		--BayesRatioFile br.tmp >>kelvin.out 2>&1
	fgrep -v Version br.tmp | perl -pe 's/Seq: 1/Seq: 3/' >> br.split
	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.split br.out .000001 || log "wordDiff for target $@ in $(PWD) FAILED"

test-USE_DL :
	echo Compilation not appropriate for this test.

test-FIXED :
	echo Fixed-grid comparison not available

clean :
	- rm -f br.out br.tmp br.split ppl.out kelvin.out*
