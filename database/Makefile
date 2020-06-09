# Kelvin db
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

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
