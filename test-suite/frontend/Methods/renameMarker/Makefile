test :
	- rm -rf case*.{locus,map,freq}
	perl -I$(KELVIN_ROOT) ./test.pl
	diff case1.freq.baseline case1.freq
	diff case2.freq.baseline case2.freq
	diff case2.map.baseline case2.map
	diff case3.locus.baseline case3.locus
	diff case3.map.baseline case3.map
	diff case4.freq.baseline case4.freq
	diff case4.locus.baseline case4.locus
	diff case4.map.baseline case4.map
	diff case5.locus.baseline case5.locus
	diff case6.locus.baseline case6.locus
	diff case6.map.baseline case6.map
	diff case7.freq.baseline case7.freq
	diff case7.locus.baseline case7.locus
	diff case7.map.baseline case7.map


clean :
	- rm -rf case*.{locus,map,freq}
