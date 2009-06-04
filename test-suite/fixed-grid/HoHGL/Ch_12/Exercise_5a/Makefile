test :
	- rm -rf *.out
	(time $(TEST_KELVIN) kelvin.EX3) >kelvin-EX3.out 2>&1
	diff -i br.out br-EX3.out-baseline
	diff -i ppl.out ppl-EX3.out-baseline
	(time $(TEST_KELVIN) kelvin.EX5a) >kelvin-EX5a.out 2>&1
	diff -i br.out br-EX5a.out-baseline
	diff -i ppl.out ppl-EX5a.out-baseline
	(time $(TEST_KELVIN) kelvin.EX5b) >kelvin-EX5b.out 2>&1
	diff -i br.out br-EX5b.out-baseline
	diff -i ppl.out ppl-EX5b.out-baseline
	(time $(TEST_KELVIN) kelvin.EX5c) >kelvin-EX5c.out 2>&1
	diff -i br.out br-EX5c.out-baseline
	diff -i ppl.out ppl-EX5c.out-baseline

clean :
	- rm -f *.out

