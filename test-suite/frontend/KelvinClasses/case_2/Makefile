test :
	- rm -rf *.new*
	perl -I $(KELVIN_ROOT) classtest.pl
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.post.new.1 ped.post-test1 0
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.post.new.2 ped.post-test2 0
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.pre.new.1 ped.pre-test1 0
	# This sorting bit is temporary
	cat ped.pre.new.2 | sort -k1,1g -k2,2g >ped.pre.new.2.sorted
	perl $(KELVIN_ROOT)/utils/wordDiff.pl ped.pre.new.2.sorted ped.pre-test2 0

clean :
	- rm -rf *.new*
