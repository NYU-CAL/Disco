MAKEFILE_H5 = $(PWD)/../Makefile_dir.in
include $(MAKEFILE_H5)

EXE  = vdisco

GLUT = /usr/X11R6
HDF5 = $(H55)

LIBS = 
CFLAGS = 

#LIBS = -L$(GLUT)/lib -L$(HDF5)/lib -framework GLUT -framework OpenGL -lXmu -lXext -lXi -lX11 -lm -lhdf5 -lz
UNAME = $(shell uname -s)
ifeq ($(UNAME), Linux)
	LIBS = -L$(HDF5)/lib -L$(GLUT)/lib -lhdf5 -lm
	CFLAGS += -D LINUX
endif
ifeq ($(UNAME), Darwin)
	LIBS = -L$(HDF5)/lib -framework GLUT -framework OpenGL -lhdf5 -lm
	CFLAGS += -D OSX
endif

INCL = -I$(GLUT)/include -I$(HDF5)/include
LIBS += -lhdf5 -lm

CC = g++

default: $(EXE)

$(EXE): $(EXE).c
	$(CC) -c colorbar.c
	$(CC) $(CFLAGS) -O3 -Wall -o $@ $< colorbar.o $(INCL) $(LIBS)

clean:
	rm -f *.o $(EXE)
