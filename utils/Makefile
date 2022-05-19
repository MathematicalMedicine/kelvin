# Kelvin utils
# Copyright (C) 2006, 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

HASH_INCS = hashtab.h lookupa.h recycle.h standard.h
HASH_OBJS = hashtab.o lookupa.o recycle.o
NOGSL_INCS = dists.h cdflib.h
NOGSL_OBJS = dcdflib.o dists.o ipmpar.o

INCS = utils.h sw.h polynomial.h polynomial_internal.h tpl.h sSDHandler.h pageManagement.h $(HASH_INCS) $(NOGSL_INCS)
OBJS = utils.o sw.o polynomial.o tpl.o sSDHandler.o pageManagement.o $(HASH_OBJS) $(NOGSL_OBJS)

all : ../lib/libklvnutls.a # ../lib/libklvnutls.so # pt sw sSDHandler

pt : pt.c polynomial.c polynomial.h polynomial_internal.h sw.c dists.c $(NOGSL_INCS) $(NOGSL_OBJS) $(HASH_INCS) $(HASH_OBJS)
	$(CC) pt.c -o pt polynomial.c sw.c $(NOGSL_OBJS) $(HASH_OBJS) $(CFLAGS) -DPOLYUSE_DL -DPOLYCODE_DL -ldl -lm -lpthread

sw : sw.c sw.h
	$(CC) -DMAIN sw.c -o sw $(CFLAGS) -lpthread -lm

sSDHandler : sSDHandler.c sSDHandler.h
	$(CC) -DMAIN sSDHandler.c -o sSDHandler

listener : listener.c sw.h
	$(CC) listener.c -o listener -lmysqlclient -I/usr/local/mysql/include -L/usr/local/mysql/lib \
	-I/usr/include/mysql -L/usr/lib64/mysql

../lib/libklvnutls.a : libklvnutls.a
	cp libklvnutls.a ../lib

../lib/libklvnutls.so : libklvnutls.so
	cp libklvnutls.so ../lib

libklvnutls.a : $(OBJS)
	ar -rsc libklvnutls.a $(OBJS)

libklvnutls.so : $(OBJS)
	$(CC) -shared -o libklvnutls.so $(OBJS) $(CFLAGS) $(LDFLAGS)

# The -fno-common bit is just the old PowerMac way of saying PIC
$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c -fno-common -fPIC $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -f pt sw listener sSDHandler $(OBJS) libklvnutls.a ../lib/libklvnutls.a libklvnutls.so ../lib/libklvnutls.so
