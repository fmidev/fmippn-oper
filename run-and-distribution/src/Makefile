COMPILER = gcc
CFLAGS =  -O2 -std=c99 -Wall -DLINUX -MMD -D_GNU_SOURCE -D_REENTRANT
LIBDIRS = -L. -L./lib_local -L/var/opt/lib
INCDIRS = -I./include_local -I/var/opt/include 
LIBS = -lhdf5_hl -lhdf5 -lz -lsz -lm
GDLIBS = -lm -lgd
GDALLIBS = -lhdf5_hl -lhdf5 -lz -lsz -lm -lproj -lgdal
PTHLIBS = -lpthread -lhdf5_hl -lhdf5 -lz -lsz -lm 
LDOPTS =

CC = $(COMPILER) $(INCDIRS) $(CFLAGS)
EXES = acc1h

all	: $(EXES)

acc1h:	acc1h.o
	$(CC) acc1h.o -o acc1h $(LDOPTS) -lm

replace_grey:	replace_grey.o
	$(CC) replace_grey.o -o replace_grey $(LDOPTS)

install	:
		@\cp $(EXES) ~/ppn/fmippn-postprocess/bin
clean	:
		@\rm -f *.o *.d
