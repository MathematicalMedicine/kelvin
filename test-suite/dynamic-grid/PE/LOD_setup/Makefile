#
# There are two clients which request all the same pedigrees but different
# sets of loci. There are two servers, which provide different pedigrees
# on different maps for all loci.
#

test :
	echo "Nope"
#	- rm -rf *.out
#	(time $(TEST_KELVIN) kelvin.conf) >kelvin.out 2>&1
#	perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-baseline .000001

test-USE_DL :
	echo "Nope"

test-FIXED :
	-rm -f LRTest*
	(perl $(KELVIN_ROOT)/LRTest.pl) >kelvin-FIXED.out 2>&1

clean :
	- rm -rf *.out *.so *memory.dat compiled LRTest*

