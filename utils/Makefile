# Kelvin
# Copyright 2008, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = utils.h sw.h tpl.h hashtab.h lookupa.h recycle.h standard.h sSDHandler.h dists.h polynomial.h tools.h
OBJS = utils.o sw.o tpl.o hashtab.o lookupa.o recycle.o sSDHandler.o dcdflib.o dists.o ipmpar.o polynomial.o tools.o

all : ../lib/libutils.a

../lib/libutils.a : libutils.a
	cp libutils.a ../lib

libutils.a : $(OBJS)
	ar -rsc libutils.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -f $(OBJS) libutils.a ../lib/libutils.a
