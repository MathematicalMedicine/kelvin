# Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.
CC=gcc
INC=-I. -I../
# The D_REENTRANT bit is for getopt_r with pthreads under Solaris
CFLAGS=-g -D_REENTRANT
LDFLAGS=-lm
#LDFLAGS=-g -lm  # If you want to run in a debugger

calc_updated_ppl : calc_updated_ppl.o map.o iplist.o positionlist.o multidim.o
	$(CC) -o $@ calc_updated_ppl.o map.o iplist.o positionlist.o multidim.o $(LDFLAGS)

%.o : %.c 
	$(CC) $(INC) $(CFLAGS) -c -o $@ $<

clean :
	- rm -f *.o calc_updated_ppl
