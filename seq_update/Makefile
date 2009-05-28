CC=cc
INC=-I. -I../
CFLAGS=-g
LDFLAGS=-lm
#LDFLAGS=-g -lm  # If you want to run in a debugger

calc_updated_ppl : calc_updated_ppl.o map.o ippl.o multidim.o
	$(CC) -o $@ $(LDFLAGS) calc_updated_ppl.o map.o ippl.o multidim.o 

%.o : %.c 
	$(CC) $(INC) $(CFLAGS) -c -o $@ $<

clean :
	- rm -f *.o calc_updated_ppl
