.SUFFIXES: .o .f90

#compiler
FC = gfortran

#compile flags
FFLAGS = -O3

#object files created from source
SRCS = Quake.o

#rule to make .o from .f90
.f90.o:
        $(FC) $(FFLAGS) $*.f90 -c

#default to make all programs
all: Quake

#how to make program
Quake: $(SRCS)
        $(FC) $(FFLAGS) $^ -o $@

#clean
clean:
        rm -f *.o Quake
