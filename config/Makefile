# Kelvin
# Copyright (C) 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

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
