# Make file for map analisys (mapan) project
# Written by Val Kalatsky
# 2001-05-10
VERSION = 1
CC = g77
LIBS = -lcpgplot -lpgplot -lX11 -lpng
LIBDIRS = -L/usr/X11R6/lib -L$(HOME)/lib
INCLUDEDIRS = -I$(HOME)/include
MAINC = mapanm2$(VERSION).c
HFILES = mapanm2$(VERSION).h
CFLAGS = -O3 -Wall -m32
OBJS =  analysis$(VERSION).o average$(VERSION).o  \
	display$(VERSION).o getoptions$(VERSION).o \
	graphics$(VERSION).o initialization$(VERSION).o \
	interactive$(VERSION).o

mapanm2 : $(MAINC) $(OBJS) $(HFILES)
	g77 $(MAINC) $(OBJS) -o mapanm2$(VERSION) $(LIBS) $(LIBDIRS) $(INCLUDEDIRS) -m32

mapanm2o : $(MAINC) $(OBJS) $(HFILES)
	g77 $(MAINC) $(OBJS) -o mapanm2$(VERSION)o $(LIBS) $(LIBDIRS) $(INCLUDEDIRS) -O3 -inline -m32

analysis$(VERSION).o : analysis$(VERSION).c $(HFILES)
	$(CC) -c analysis$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

average$(VERSION).o : average$(VERSION).c $(HFILES)
	$(CC) -c average$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

display$(VERSION).o : display$(VERSION).c $(HFILES)
	$(CC) -c display$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

getoptions$(VERSION).o : getoptions$(VERSION).c $(HFILES)
	$(CC) -c getoptions$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

graphics$(VERSION).o : graphics$(VERSION).c $(HFILES)
	$(CC) -c graphics$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

initialization$(VERSION).o : initialization$(VERSION).c $(HFILES)
	$(CC) -c initialization$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)

interactive$(VERSION).o : interactive$(VERSION).c $(HFILES)
	$(CC) -c interactive$(VERSION).c $(INCLUDEDIRS) $(CFLAGS)
