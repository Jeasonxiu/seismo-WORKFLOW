#Directories
PWD    = $(shell pwd)
OBJDIR = ../../obj
SRCDIR = $(PWD)
#LIBDIR = $(PWD)/UTILS/lib
#IDIR   = $(PWD)/UTILS/include
#EVALLIB=/home/lei/ASDF/evalresp/.libs

CFLAGS = -check

# new version with default sac libraries
SACLIBDIR = /home/jas11/lib
LIBS = -lsacio -lsac
#LIB = -L/opt/seismo/lib -lDRWFiles -lf90recipes -lDSacio -lDSacLib -lSacTools -lm

############################
#compiler option
#OPT = -I${SHARED}
#OPT = -std03
FC = ifort
CC = icc
MPIFC = mpif90
MPICC = mpicc
CFLAGS= -g -O0 -check all

_OBJ = ma_constants.o ma_variables.o ascii_rw.o \
			 ma_sub2.o ma_sub.o measure_adj.o \
			 mpi_weighting_subs.o rotate_subs.o \
			 measure_adj_interface_subs.o measure_adj_interface.o

OBJ = $(patsubst %, ${OBJDIR}/%, $(_OBJ))

## set ADIOS_DIR here or before doing make
#override ADIOS_DIR:=/home/jas11/adios-1.5.0
#override ADIOS_INC:=` ${ADIOS_DIR}/bin/adios_config -c -f`
#override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`
ADIOS_INC=$(shell adios_config -cf)
ADIOS_FLIB=$(shell adios_config -lf)

##########################################################
PROG = TEST
default: ${PROG}

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)
.F.o:
	$(OBJDIR)/%.o: $(SRCDIR)/%.F
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

#$(OBJDIR)/%.o: $(SRCDIR)/%.c
#	  $(MPICC) -c -o $@ $< 

${PROG}: $(OBJ)


.PHONY:clean print_var cleanall

print_var:
	@echo $(OBJ)
	@echo $(SRCDIR)

clean:
	rm -f  ${LIB_ALL} ${PROGS} *.o *.mod *.a $(OBJDIR)/*

cleanall:
	rm -f  iasp91.*
	cd ${TAULIBDUR} ; make -f make_gfortran clean

