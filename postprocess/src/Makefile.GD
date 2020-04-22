COMPILER=gcc
CFLAGS=-O -DHAVE_XPM -Wall
GDLIBS=-lm -lgd
# GDLIBS=-lm -lgd -lpng -lz -lX11
MATHLIBS=-lm
INCLUDEDIRS=-I./gd -I/usr/include/gd -I/usr/include/X11 -I/usr/X11R6/include/X11
LIBDIRS=-L. -L/usr/lib64

CC=$(COMPILER) $(INCLUDEDIRS)
LINK=$(CC) $(LIBDIRS) $(LIBS)
  
GD_SRCS = 	gdmap.c

GD_OBJS =	${GD_SRCS:.c=.o}

EXES =		${GD_SRCS.c=}

all:		${EXES}

clean_all:;	rm -f ${EXES} *.o

gdmap:		gdmap.o
		$(CC) gdmap.o -o gdmap $(LIBDIRS) $(GDLIBS)  
