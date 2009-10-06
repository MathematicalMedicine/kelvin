# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = utils.h sw.h polynomial.h polynomial_internal.h tpl.h hashtab.h lookupa.h recycle.h standard.h sSDHandler.h dists.h pageManagement.h
OBJS = utils.o sw.o polynomial.o tpl.o hashtab.o lookupa.o recycle.o sSDHandler.o dcdflib.o dists.o ipmpar.o pageManagement.o

all : ../lib/libutils.a ../lib/libklvnutls.so # pt sw sSDHandler

pt : pt.c polynomial.c polynomial.h polynomial_internal.h
	gcc pt.c -o pt polynomial.c sw.c dists.c $(CFLAGS) -DUSE_GSL -DPOLYUSE_DL -DPOLYCODE_DL -ldl -lgsl -lgslcblas -lm

sw : sw.c sw.h
	gcc -DMAIN sw.c -o sw $(CFLAGS) -lm

sSDHandler : sSDHandler.c sSDHandler.h
	gcc -DMAIN sSDHandler.c -o sSDHandler

listener : listener.c sw.h
	gcc listener.c -o listener -lmysqlclient -I/usr/local/mysql/include -L/usr/local/mysql/lib \
	-I/usr/include/mysql -L/usr/lib64/mysql

../lib/libutils.a : libutils.a
	cp libutils.a ../lib

../lib/libklvnutls.so : libklvnutls.so
	cp libklvnutls.so ../lib

libutils.a : $(OBJS)
	ar -rsc libutils.a $(OBJS)

libklvnutls.so : $(OBJS)
	gcc -shared -o libklvnutls.so $(OBJS) $(CFLAGS) $(LDFLAGS)

$(OBJS) : %.o : %.c $(INCS)
	# The -fno-common bit is just the old PowerMac way of saying PIC
	$(CC) -c -fno-common -fPIC $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -f pt $(OBJS) libutils.a ../lib/libutils.a libklvnutls.so ../lib/libklvnutls.so
