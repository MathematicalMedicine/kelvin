# Kelvin
# Copyright 2010, The Research Institute At Nationwide Children's Hospital
# Permission granted to distribute and use for non-profit educational purposes only.

# Invoked by main Kelvin Makefile. You'll get undefined symbols otherwise.

INCS = config.h model_options.h model_type.h model_range.h lambda_cell.h

OBJS = config.o read_pedfile.o model_range.o model_type.o

all : ../lib/libconfig.a

clean :
	rm -f $(OBJS) libconfig.a ../lib/libconfig.a

../lib/libconfig.a : libconfig.a
	cp libconfig.a ../lib

libconfig.a : $(OBJS)
	ar -rsc libconfig.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@
