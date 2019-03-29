# Lets make the RADSMM library
#CC	= cc
#CFLGS  = -O -Wall

OBJS  = ./RADSMM_setup.o \
	./RADSMM_index.o \
	./RADSMM_seek.o \
	./RADSMM_read.o \
	./RADSMM_write.o \
	./RADSMM_file.o \
	./RADSMM_create_file.o \
	./RADSMM_label.o

######################################################################
# Default target
.PHONY: install
install: libRADSMM.a
	-@cp {RADSMM,RADSMM_file_header}.h $(INCDIR)
	-@cp libRADSMM.a $(LIBDIR)

######################################################################
libRADSMM.a:	$(OBJS)
	-@rm -f libRADSMM.a
	ar -rc libRADSMM.a $(OBJS)

clean:
	-@rm -f *.o
	-@rm -f *.a
# some rules
.c.o:
	$(CC) $(CFLGS) -c $< -o $@

.SUFFIXES:
.SUFFIXES:	.o .c .f

# dependences

RADSMM_setup.o :	RADSMM_setup.c RADSMM.h
RADSMM_index.o :	RADSMM_index.c RADSMM.h
RADSMM_create_file.o :	RADSMM_create_file.c RADSMM.h RADSMM_file_header.h
RADSMM_read.o :		RADSMM_read.c RADSMM.h
RADSMM_write.o :	RADSMM_write.c RADSMM.h
RADSMM_seek.o :		RADSMM_seek.c RADSMM.h
RADSMM_file.o :		RADSMM_file.c RADSMM.h RADSMM_file_header.h
RADSMM_label.o :	RADSMM_label.c RADSMM.h
