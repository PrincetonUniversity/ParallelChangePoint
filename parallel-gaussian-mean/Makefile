# Makefile for parallel gaussian mean changepoint detection
#
# (C) Haw Yang 2012
# 

UNAME := $(shell uname)
PROGRAM = Parallel-GM
SOURCES = CPTree.c C_ac.c CheckCP.c FindCP.c MakeCPArray.c main.c \
          Vostrikova.c zbrent.c changepoint.h critical_values.h
OBJECTS = CPTree.o C_ac.o CheckCP.o FindCP.o MakeCPArray.o main.o \
			Vostrikova.o zbrent.o
LIB = -lm -lgsl -lgslcblas -lrt
CFLAGS = -O3
CC = mpicc

MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

ifeq ($(strip $(OSTYPE)),msys)
  ########## MinGW - msys ##########
  # For Windows, make sure to include the libgsl.dll path in the system PATH
  # variable. This can be done by right-clicking 'My Computer' under 
  # 'Properties.'
  PLATFORM = MinGW
  INCDIRFLAG = -I/c/Program\ Files/GnuWin32/include
  LIBDIRFLAG = -L/c/Program\ Files/GnuWin32/lib
endif
ifeq ($(UNAME), Darwin)
  ########## Mac OS X ##########
  PLATFORM = Darwin
  INCDIRFLAG = -I/opt/local/include
  LIBDIRFLAG = -L/opt/local/lib
endif
ifeq ($(strip $(OSTYPE)),linux-gnu)
  ########## Linux - linux-gnu ##########
  PLATFORM = Linux
endif

$(PROGRAM): $(OBJECTS)
	$(CC) $(OBJECTS) $(MPI_LINK_FLAGS) -o $(PROGRAM) $(LIB) $(LIBDIRFLAG) 

.c.o:   ; $(CC) $(MPI_COMPILE_FLAGS) -c $*.c $(INCDIRFLAG)

clean:
	rm -f $(PROGRAM) $(OBJECTS)
