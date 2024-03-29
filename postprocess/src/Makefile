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
EXES = ensmean_determweighted thread_member_interp prob_thresholding

all	: $(EXES)

ensmean:	ensmean.o
		$(CC) ensmean.o -o ensmean $(LDOPTS) $(LIBDIRS) $(LIBS)

gdmap:	gdmap.o
		$(CC) gdmap.o -o gdmap $(LDOPTS) $(LIBDIRS) $(GDLIBS)

ensmean_determweighted:	ensmean_determweighted.o
		$(CC) ensmean_determweighted.o -o ensmean_determweighted $(LDOPTS) $(LIBDIRS) $(LIBS)

member_interp:	member_interp.o
		$(CC) member_interp.o -o member_interp $(LDOPTS) $(LIBDIRS) $(LIBS)

thread_member_interp:	thread_member_interp.o
		$(CC) thread_member_interp.o -o thread_member_interp $(LDOPTS) $(LIBDIRS) $(PTHLIBS)

getmotion:	getmotion.o
		$(CC) getmotion.o -o getmotion $(LDOPTS) $(LIBDIRS) $(LIBS)

prob_thresholding:	prob_thresholding.o
		$(CC) prob_thresholding.o -o prob_thresholding $(LDOPTS) $(LIBDIRS) $(GDALLIBS) 

reprojection_radardata:	reprojection_radardata.o
		$(CC) reprojection_radardata.o -o reprojection_radardata $(LDOPTS) $(LIBDIRS) $(GDALLIBS)

extract_motion:	extract_motion.o
		$(CC) extract_motion.o -o extract_motion $(LDOPTS) $(LIBDIRS) $(LIBS)

hourmask:	hourmask.o
		$(CC) hourmask.o -o hourmask $(LDOPTS)

install	:
		@\cp $(EXES) ~/ppn/fmippn-postprocess/bin
clean	:
		@\rm -f *.o *.d
