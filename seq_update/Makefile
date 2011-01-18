CC=gcc
INC=-I. -I../
# The D_REENTRANT bit is for getopt_r with pthreads under Solaris
CFLAGS=-g -D_REENTRANT
LDFLAGS=-lm
#LDFLAGS=-g -lm  # If you want to run in a debugger

calc_updated_ppl : calc_updated_ppl.o map.o iplist.o positionlist.o multidim.o
	$(CC) -o $@ $(LDFLAGS) calc_updated_ppl.o map.o iplist.o positionlist.o multidim.o 

%.o : %.c 
	$(CC) $(INC) $(CFLAGS) -c -o $@ $<

clean :
	- rm -f *.o calc_updated_ppl
