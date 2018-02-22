test :
	- rm -rf *.out
	@echo -e "\n\
    **** THIS TEST IS CURRENTLY DISABLED ****\n\
"
	#- rm -f ped.post.counted ped.post.counts kelvin.freq *.out
	# Generate a new frequency file and count files for comparison
	#(time $(KELVIN_SCRIPT) kelvin.conf --SkipAnalysis) >kelvin.out 2>&1
	#perl $(KELVIN_ROOT)/utils/wordDiff.pl kelvin.freq kelvin.freq-baseline 0
	#diff -ib ped.post.counted ped.post.counted-baseline
	#diff -ib ped.post.counts ped.post.counts-baseline
	#- rm -f ped.post.counted ped.post.counts kelvin.freq *.out
	# Now run the analysis with Kelvin and preserve the results
	#(time $(KELVIN_SCRIPT) kelvin.conf) >kelvin.out 2>&1
	#cp br.out br.out-by-Kelvin
	#cp ppl.out ppl.out-by-Kelvin
	# Finally run without the front-end and compare the results with the front-end
	#(time $(TEST_KELVIN) kelvin.conf --FrequencyFile kelvin.freq) >kelvin.out 2>&1
	#perl $(KELVIN_ROOT)/utils/wordDiff.pl br.out br.out-by-Kelvin .000001
	#perl $(KELVIN_ROOT)/utils/wordDiff.pl ppl.out ppl.out-by-Kelvin .000001

clean :
	- rm -f ped.post.counted ped.post.counts kelvin.freq* *.out

