# Make file for image analisys (iman) project
# Written by Val Kalatsky
# 2001-03-23
VERSION = 2
CC = g77
LIBS =  -lcpgplot -lpgplot -lX11 -lpng  -lbesselI -lfft
LIBDIRS = -L/usr/X11R6/lib -L/usr/local/lib
INCLUDEDIRS = -I/usr/local/include
MAINC = iman$(VERSION).c
HFILES = iman1.h fileheader.h
CFLAGS = -O3 -Wall -m32
OBJS =  getoptions$(VERSION).o initializeroutines$(VERSION).o  \
	frameroutines$(VERSION).o printroutines$(VERSION).o      \
	statistics$(VERSION).o graphics$(VERSION).o 	       \
	cocktailaverageroutines$(VERSION).o displayroutines$(VERSION).o \
	convertroutines$(VERSION).o traceanalysis$(VERSION).o \
	greenroutines$(VERSION).o synchronization$(VERSION).o \
	inverse$(VERSION).o interframetime$(VERSION).o \
	fitroutines$(VERSION).o compressor$(VERSION).o binning$(VERSION).o \
	fileheader.o

iman$(VERSION) : $(MAINC) $(OBJS) $(HFILES)
	$(CC) $(MAINC) $(OBJS) -o iman$(VERSION) $(LIBDIRS) $(INCLUDEDIRS) $(LIBS) -m32

imano : $(MAINC) $(OBJS) $(HFILES)
	g77 $(MAINC) $(OBJS) -o iman$(VERSION)o $(LIBS) $(LIBDIRS) $(INCLUDEDIRS) -O3 -m32

getoptions$(VERSION).o : getoptions$(VERSION).c $(HFILES)
	$(CC) -c getoptions$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

initializeroutines$(VERSION).o : initializeroutines$(VERSION).c $(HFILES)
	$(CC) -c initializeroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

frameroutines$(VERSION).o : frameroutines$(VERSION).c $(HFILES)
	$(CC) -c frameroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

printroutines$(VERSION).o : printroutines$(VERSION).c $(HFILES)
	$(CC) -c printroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

statistics$(VERSION).o : statistics$(VERSION).c $(HFILES)
	$(CC) -c statistics$(VERSION).c $(INCLUDEDIRS) -I/usr/local/src/functions $(CFLAGS)

graphics$(VERSION).o : graphics$(VERSION).c $(HFILES)
	$(CC) -c graphics$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

cocktailaverageroutines$(VERSION).o : cocktailaverageroutines$(VERSION).c $(HFILES)
	$(CC) -c cocktailaverageroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

displayroutines$(VERSION).o : displayroutines$(VERSION).c $(HFILES)
	$(CC) -c displayroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

convertroutines$(VERSION).o : convertroutines$(VERSION).c $(HFILES)
	$(CC) -c convertroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

traceanalysis$(VERSION).o : traceanalysis$(VERSION).c $(HFILES)
	$(CC) -c traceanalysis$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

greenroutines$(VERSION).o : greenroutines$(VERSION).c $(HFILES)
	$(CC) -c greenroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

synchronization$(VERSION).o : synchronization$(VERSION).c $(HFILES)
	$(CC) -c synchronization$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

inverse$(VERSION).o : inverse$(VERSION).c $(HFILES)
	$(CC) -c inverse$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

interframetime$(VERSION).o : interframetime$(VERSION).c $(HFILES)
	$(CC) -c interframetime$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

fitroutines$(VERSION).o : fitroutines$(VERSION).c $(HFILES)
	$(CC) -c fitroutines$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

compressor$(VERSION).o : compressor$(VERSION).c $(HFILES)
	$(CC) -c compressor$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

binning$(VERSION).o : binning$(VERSION).c $(HFILES)
	$(CC) -c binning$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

fileheader.o : fileheader.c $(HFILES)
	$(CC) -c fileheader.c $(INCLUDEDIRS) $(CFLAGS)


imagenerator : imagenerator$(VERSION).c $(HFILES)
	$(CC) -o $@ imagenerator$(VERSION).c $(CFLAGS) $(INCLUDEDIRS) $(LIBDIRS) $(LIBS)  

