# Kelvin db
# Copyright (C) 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = StudyDB.h databaseSupport.h
OBJS = databaseSupport.o

all : ../lib/libklvndb.a

../lib/libklvndb.a : libklvndb.a
	cp libklvndb.a ../lib

../lib/libklvndb.so : libklvndb.so
	cp libklvndb.so ../lib

libklvndb.a : $(OBJS)
	ar -rsc libklvndb.a $(OBJS)

libklvndb.so : $(OBJS)
	gcc -shared -o libklvndb.so $(OBJS) $(CFLAGS) $(LDFLAGS)

# The -fno-common bit is just the old PowerMac way of saying PIC
$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c -fno-common -fPIC $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -f pt $(OBJS) libklvndb.a ../lib/libklvndb.a libklvndb.so ../lib/libklvndb.so
