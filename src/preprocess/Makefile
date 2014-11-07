#Directories
PWD    = $(shell pwd)
OBJDIR = ../../obj
SRCDIR = $(PWD)
#LIBDIR = $(PWD)/UTILS/lib
#IDIR   = $(PWD)/UTILS/include
#EVALLIB=/home/lei/ASDF/evalresp/.libs

CFLAGS = -check noarg_temp_created -g -O0

# new version with default sac libraries
SACLIBDIR = /home/jas11/sac-101.6a/build/src
LIBS = -lsacio -lsac
 #LIB = -L/opt/seismo/lib -lDRWFiles -lf90recipes -lDSacio -lDSacLib -lSacTools -lm


############################
#compiler option
#OPT = -I${SHARED}
#OPT = -std03
FC = ifort 
CC = cc
MPIFC = mpif90
MPICC = mpicc
CFLAGS= -g

_OBJ = 	process_par.o process_var.o process_subs.o process_subs2.o\
				process_asdf.o

OBJ = $(patsubst %, ${OBJDIR}/%, $(_OBJ))

## set ADIOS_DIR here or before doing make
#override ADIOS_DIR:=/home/jas11/adios-1.5.0
#override ADIOS_INC:=` ${ADIOS_DIR}/bin/adios_config -c -f`
#override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`

ASDFLIBDIR=$(ASDFHOME)/lib
ASDFINCDIR=$(ASDFHOME)/include

ADIOS_INC=$(shell adios_config -cf)
ADIOS_FLIB=$(shell adios_config -lf)

##########################################################
PROG = TEST
default: ${PROG}

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC) -I$(ASDFINCDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC) -I$(ASDFINCDIR)
.F.o:
	$(OBJDIR)/%.o: $(SRCDIR)/%.F
	  $(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC) -I$(ASDFINCDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	  $(MPICC) -c -o $@ $< 

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

