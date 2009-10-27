test :
	- rm -rf PC* BR.Out-PCs ppl.out kelvin.out*
	$(TEST_KELVIN) kelvin.conf >kelvin.out 2>&1
	$(KELVIN_ROOT)/PedCount.pl -config kelvin.conf -split=1 -subdir -write
	$(TEST_KELVIN) kelvin.conf-PC1 >kelvin.out-PC1 2>&1
	cat PC1/BR.Out >BR.Out-PCs
	$(TEST_KELVIN) kelvin.conf-PC2 >kelvin.out-PC2 2>&1
	perl -p -i -e 's/Seq: 1/Seq: 2/' PC2/BR.Out
	grep -v Version PC2/BR.Out >>BR.Out-PCs
	$(TEST_KELVIN) kelvin.conf-PC3 >kelvin.out-PC3 2>&1
	perl -p -i -e 's/Seq: 1/Seq: 3/' PC3/BR.Out
	grep -v Version PC3/BR.Out >>BR.Out-PCs
	diff -i BR.Out-PCs br.out

test-USE_DL :
	echo Compilation not appropriate for this test.

test-FIXED :
	echo Fixed-grid comparison not available

clean :
	- rm -rf PC* BR.Out-PCs ppl.out kelvin.out* LRTest*
