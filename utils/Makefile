# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = utils.h sw.h tpl.h hashtab.h lookupa.h recycle.h standard.h sSDHandler.h dists.h pageManagement.h
OBJS = utils.o sw.o tpl.o hashtab.o lookupa.o recycle.o sSDHandler.o dcdflib.o dists.o ipmpar.o pageManagement.o

all : ../lib/libutils.a pt sw sSDHandler

polynomial.o : polynomial.c polynomial.h
	gcc polynomial.c -c $(CFLAGS) $(INCFLAGS) -o polynomial.o

pt : pt.c polynomial.o sw.o
	gcc polynomial.c -c -o polynomial.o -DPOLYUSE_DL -DPOLYCODE_DL
	gcc pt.c -o pt polynomial.o -L../lib -lutils -lm -ldl
	rm polynomial.o

sw : sw.c sw.h
	gcc -DMAIN sw.c -o sw -lm

sSDHandler : sSDHandler.c sSDHandler.h
	gcc -DMAIN sSDHandler.c -o sSDHandler

listener : listener.c sw.h
	gcc listener.c -o listener -lmysqlclient -I/usr/local/mysql/include -L/usr/local/mysql/lib \
	-I/usr/include/mysql -L/usr/lib64/mysql

../lib/libutils.a : libutils.a
	cp libutils.a ../lib

libutils.a : $(OBJS) polynomial.o
	ar -rsc libutils.a $(OBJS) polynomial.o

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -f pt $(OBJS) polynomial.o libutils.a ../lib/libutils.a
