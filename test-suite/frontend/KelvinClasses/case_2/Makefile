test :
	- rm -rf *.new*
	perl -I $(KELVIN_ROOT) classtest.pl
	diff -bB ped.post.new.1 ped.post-test1
	diff -bB ped.post.new.2 ped.post-test2
	diff -bB ped.pre.new.1 ped.pre-test1
	diff -bB ped.pre.new.2 ped.pre-test2

clean :
	- rm -rf *.new*
