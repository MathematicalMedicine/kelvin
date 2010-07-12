test :
	- rm -f ped.post.counted ped.post.counts *.out
	(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	diff -ib ped.post.counted ped.post.counted-baseline
	diff -ib ped.post.counts ped.post.counts-baseline

clean :
	- rm -f ped.post.counted ped.post.counts *.out

