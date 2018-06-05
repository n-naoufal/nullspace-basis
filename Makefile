include make.inc

#######################################################################
#  This makefile creates the Fortran example interface to use the
#  C routines in SuperLU.
#######################################################################

HEADER   = ../SRC
LIBS	= $(SUPERLULIB) $(BLASLIB) -lm

# double real
NULLBASIS	= nullspace_basis.o dsorti.o c_fortran_dgssv.o

all:	nullbasis 

nullbasis: $(NULLBASIS) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(nullbasis) $(LIBS) -o $@

nullspace_basis.o: nullspace_basis.f90
	$(FORTRAN) $(LOADOPTS) -c $< 

dsorti.o: dsorti.f90
	$(FORTRAN) $(LOADOPTS) $(LIBS) -c $<

c_fortran_dgssv.o: c_fortran_dgssv.c
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o nullbasis 

