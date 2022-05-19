# Kelvin
# Copyright (C) 2006, 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

INCS = pedigree.h locus.h genotype_elimination.h \
       likelihood.h pedlib.h allele_set_recoding.h

OBJS = locus.o allele_set_recoding.o \
       genotype_elimination.o parental_pair.o likelihood.o

all : ../lib/libped.a

clean :
	rm -f $(OBJS) libped.a ../lib/libped.a

../lib/libped.a : libped.a
	cp libped.a ../lib

libped.a : $(OBJS)
	ar -rsc libped.a $(OBJS)

$(OBJS) : %.o : %.c $(INCS)
	$(CC) -c $(CFLAGS) $(INCFLAGS) $< -o $@
