#**************************************************************
#   Copyright (C) 2022
#   FastOpt GmbH, Hamburg, Germany
#   Email: info@FastOpt.de
#   All rights reserved.
#
#**************************************************************

RM		= rm -rf
#does not work with tcsh#   SHELL += -o pipefail
RUN = nohup # could be: mpirun 

#-----------------------------------------------------------------
# TAF settings
#-----------------------------------------------------------------
TAF		= /scratch/ifmrs/u24331/Important/staf -version 5.9.7  #staf
TAFFLAGS	= -input x -output fc -toplevel model
#test# TAFFLAGS	= -input x -output fc -toplevel master
#TAFFLAGS	+= -keep
#TAFFLAGS	+= -pure
#TAFFLAGS	+= -pureadd
TAFFLAGS	+= -f08 -arglist
#TAFFLAGS        += -warn_step
#TAFFLAGS        += -mpi
#TAFFLAGS	+= -newv
#TAFFLAGS       += -showrecomp

#-----------------------------------------------------------------
# compiler settings
#-----------------------------------------------------------------
#FC = gfortran
FC =ifx
ifneq (,$(findstring gfortran,$(FC))) # if search for "gfortran" in  FC is not empty
  DOUBLE  = -fdefault-real-8 -fdefault-double-8
  FFLAGS  = -J $(BLDDIR) -g
  OPTFLAGS = -O2
#DBGFLAGS += -g -ftrapuv -traceback -init=snan -init=arrays -fpe0 # debugging options
  DBGFLAGS  = -g -Og -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow -fsanitize=bounds-strict -fbacktrace -Wall
else ifneq (,$(findstring ifx,$(FC))) # if search for "ifort" in  FC is not empty
  DOUBLE  = -r8
  FFLAGS  = -module $(BLDDIR)
  OPTFLAGS = -O -g
  DBGFLAGS = -g -ftrapuv -traceback -init=snan -init=arrays -fpe0 -check all  # debugging options
#  DBGFLAGS += -gen-interfaces -warn interfaces -warn shape
endif

FFLAGS += $(DOUBLE) # compile in double precision
TAFFLAGS        += -r8 # assume default reals are double precision

ifdef DEBUG # if "DEBUG" is set, use DBGFLAGS flags, otherwise use OPTFLAGS:
  ifeq "$(DEBUG)" "0"
    FFLAGS        += $(OPTFLAGS)
  else	
    FFLAGS        += $(DBGFLAGS)
  endif
  ifeq "$(DEBUG)" "2" # DEBUG=2 sets macro to allow for code variants, usable with .F90 files
    FCPPFLAGS+=-DDEBUG
  endif
else # default:
  FFLAGS        += $(OPTFLAGS)
endif
# for pumax_stub.c
CC = gcc
CFLAGS = -O3 -I /usr/lib/X11/include

#-----------------------------------------------------------------
# paths
#-----------------------------------------------------------------
DRVDIR = ./drivers
SRCDIR = ./src
BLDDIR = ./bld
RUNDIR = ./run
TAPEDIR = $(RUNDIR)/taftapes
BINDIR=./bin

$(BLDDIR) $(RUNDIR) $(TAPEDIR):
	mkdir -p $@

PLASIMSRC = resmod.F90        \
            plasimmod.F90     \
            calmod.F90        \
			oceanmod.F90      \
			miscmod.F90       \
			landmod.F90       \
            simba.F90         \
			surfmod.F90       \
            radmod.F90        \
			srvio.F90 		  \
			observation.F90	  \
			tafmod.F90        \
            icemod.F90        \
            seamod.F90        \
            rainmod.F90       \
            fluxmod.F90       \
            restartmod.F90    \
            wrapper.F90       \
            mpimod_stub.F90   \
            fftmod_flow_directives.F90  \
            fftmod.F90        \
            cpl_stub.F90      \
            guimod_stub.F90   \
            p_earth.F90       \
            plasim.F90        \
            gaussmod.F90      \
            legmod.F90        \
            outmod.F90        \
            lsgmod.F90		  \
			
SRC = $(PLASIMSRC) # add further source files of $(SRCDIR) here in order of dependency

#additional:
SRCC =  pumax_stub.c

#-----------------------------------------------------------------
# define names after preprocessing with cpp and TAF
#-----------------------------------------------------------------
F90SRC	= $(SRC:.F90=.f90) # names of preprocessed sources
BLDSRC	= $(addprefix $(BLDDIR)/,$(SRC:.F90=.f90)) # names with path
BLDOBJ  = $(BLDSRC:.f90=.o) # names of object files
BLDOBC	= $(addprefix $(BLDDIR)/,$(SRCC:.c=.o)) # names of object files with source in c

TLMSRC  = $(BLDSRC:.f90=_tl.f90) # names of files generated for TLM, with path
TLMOBJ  = $(BLDSRC:.f90=_tl.o)
ADMSRC  = $(BLDSRC:.f90=_ad.f90) $(BLDDIR)/taf_ad.f90  # names of files generated for ADM, with path; may need to add taf_ad.f90
ADMOBJ  = $(BLDSRC:.f90=_ad.o) $(BLDDIR)/taf_ad.o
DOTSRC  = $(BLDSRC:.f90=_dot.f90) $(BLDDIR)/taf_dot.f90  # names of files generated for DOT, with path; may need to add taf_dot.f90
DOTOBJ  = $(BLDSRC:.f90=_dot.o) $(BLDDIR)/taf_dot.o

#-----------------------------------------------------------------
# show valid targets
#-----------------------------------------------------------------
help:
	@echo 'Targets: (all prerequisites are built automatically)'
	@echo 'make cost 	: evaluate function'
	@echo 'make cost2	: evaluate function twice'
	@echo 'make tlm 	: generate tangent linear'
	@echo '` 	: generate adjoint'
	@echo 'make checktlm	: check tangent'
	@echo 'make checkadm	: check adjoint'
	@echo 'make run{cost|cost2|checktlm|checkadm}: run the executable'

.phony: help tlm adm cost cost2 checktlm checkadm clean cleanrun adjoint

#-----------------------------------------------------------------
# generate derivative code
#-----------------------------------------------------------------
$(TLMSRC): taf_tlm.log
	:
$(ADMSRC): taf_adm.log
	:
$(DOTSRC): taf_dot.log
	:

nrm: $(BLDSRC) | $(BLDDIR) #$(BLDDIR)/flow_tl.f90
	cd $(BLDDIR) && $(TAF) -l - $(TAFFLAGS) $(TAFDIFF) -lint $(F90SRC) > $(PWD)/taf_nrm.log 2>&1



tlm taf_tlm.log: $(BLDSRC) | $(BLDDIR) #$(BLDDIR)/flow_tl.f90
#not valid for sh shell: 	cd $(BLDDIR) && $(TAF) -l - $(TAFFLAGS) $(TAFDIFF) -forward $(F90SRC) >& $(PWD)/taf_tlm.log
	cd $(BLDDIR) && $(TAF) -l - $(TAFFLAGS) $(TAFDIFF) -forward $(F90SRC) > $(PWD)/taf_tlm.log 2>&1


adm taf_adm.log: $(BLDSRC) | ${BLDDIR} ${TAPEDIR} #$(BLDDIR)/flow_ad.f90
	@echo ''
	@echo 'shell=' $(SHELL)
	@echo ''
	@echo 'adm taf_adm.log'
	@echo 'BLDSRC  ' $(BLDSRC)
	@echo ''
	@echo 'BLDSRC  ' $(BLDSRC)
	@echo ''
	@echo 'BLDDIR  ' $(BLDDIR)
	@echo ''
	@echo 'TAF     ' $(TAF)
	@echo 'TAFFLAGS' $(TAFFLAGS)
	@echo 'TAFDIFF ' $(TAFDIFF)
	@echo 'F90SRC  ' $(F90SRC)
	@echo ''
	@echo 'PWD     ' $(PWD)
#not valid for sh shell: cd $(BLDDIR) && $(TAF) -l - $(TAFFLAGS) $(TAFDIFF) -reverse $(F90SRC) >& $(PWD)/taf_adm.log
	cd $(BLDDIR) && $(TAF) -l - $(TAFFLAGS) $(TAFDIFF) -reverse $(F90SRC) > $(PWD)/taf_adm.log 2>&1

dot taf_dot.log: $(BLDSRC) | ${BLDDIR} ${TAPEDIR} #$(BLDDIR)/flow_ad.f90
	cd $(BLDDIR) && $(TAF) -l - $(TAFFLAGS) $(TAFDIFF) -dottest -forward -reverse -admark _dot -dotonly "model,master,gridpointd,gridpointa,surfstep,rainstep,mkdqtgp,kuo,mkshallow,mkdca,mklsp,mkclouds,mkrain,seastep,icestep,oceanstep" $(F90SRC) > $(PWD)/taf_dot.log 2>&1

#-----------------------------------------------------------------
# generate, compile, and run code
#-----------------------------------------------------------------
cost: $(BLDOBJ) $(BLDOBC)  $(DRVDIR)/prgcost.F90 | ${BLDDIR} 
	$(LINK.F) $(FPPFLAGS) $^ -o $@
	mv cost $(BINDIR)
cost2: $(BLDOBJ) $(BLDOBC) $(DRVDIR)/prgcost2.F90 | ${BLDDIR} 
	@echo ''
	@echo 'FC' $(FC)
	@echo 'FFLAGS' $(FFLAGS)
	@echo 'BLDOBC' $(BLDOBC)
	@echo ''
	@echo 'F90SRC' $(F90SRC)
	@echo ''
	@echo '$^' $^
	@echo '$@' $@
	@echo 'LINK.F  =' $(LINK.F)
	@echo 'FPPFLAGS=' $(FPPFLAGS)
	@echo 'CC' $(CC)
	@echo 'CFLAGS' $(CFLAGS)
	@echo ''
	$(LINK.F) $(FPPFLAGS)  $^ -o $@
	mv cost2 $(BINDIR)


checktlm: $(TLMOBJ) $(BLDOBC) $(DRVDIR)/prgchecktl.F90 | ${BLDDIR} 
	$(LINK.F) $(FPPFLAGS)  $^ -o $@

checkadm: $(ADMOBJ) $(BLDOBC) $(DRVDIR)/prgcheckad.F90 | ${BLDDIR} 
	$(LINK.F) $(FPPFLAGS)  $^ -o $@
	mv checkadm $(BINDIR)

adjoint: $(ADMOBJ) $(BLDOBC) $(DRVDIR)/prgad.F90 | ${BLDDIR} 
	$(LINK.F) $(FPPFLAGS)  $^ -o $@
	mv adjoint $(BINDIR)

dottest: $(DOTOBJ) $(BLDOBC) $(DRVDIR)/prgdot.F90 | ${BLDDIR} 
	$(LINK.F) $(FPPFLAGS)  $^ -o $@


run%: % | ${RUNDIR}
	cd $(RUNDIR); [ -s $(PWD)/run_orig_2steps/check.par ] && cp $(PWD)/run_orig_2steps/check.par ./ || true
	cd $(RUNDIR); for FILE in $(PWD)/run_orig_2steps/*_namelist;do [ -s $${FILE} ] || cp -p $${FILE}; done
	cd $(RUNDIR); for FILE in $(PWD)/run_orig_2steps/*;do [ -s `basename $${FILE}` ] || ln -sf $${FILE}; done
	cd $(RUNDIR); ln -sf $(PWD)/y_mask_01_22.srv
	cd $(RUNDIR); $(RUN) $(PWD)/./$< | tee $(@).out

#-----------------------------------------------------------------
# rules for the compiler
#-----------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .f90 .F90
$(BLDDIR)/%.o: $(SRCDIR)/%.F90
	@[ ! -d $(BLDDIR) ] && [ ! -L $(BLDDIR) ] && mkdir -p $(BLDDIR) || true 
	$(FC) $(FFLAGS) $(FPPFLAGS) -c $< -o $@
$(BLDDIR)/%.f90: $(SRCDIR)/%.F90
	@[ ! -d $(BLDDIR) ] && [ ! -L $(BLDDIR) ] && mkdir -p $(BLDDIR) || true 
	cpp -P -traditional $(FPPFLAGS) $< >$@
$(BLDDIR)/%.f90: $(SRCDIR)/%.f90
	cp -p $< $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

#for pumax_stub.c:
$(BLDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

#-----------------------------------------------------------------
# clean up
#-----------------------------------------------------------------
clean: 
	-$(RM) $(BLDDIR)/*~
	-$(RM) $(BLDDIR)/*.o
	-$(RM) $(BLDDIR)/*.mod 

scratch: clean
	-$(RM) $(BLDDIR)/*.f90
	-$(RM) $(BLDDIR)/taf_command $(BLDDIR)/taf_output
	-$(RM) $(RUNDIR)/*
	-$(RM) $(BINDIR)/*

cleanrun:
	-$(RM) $(RUNDIR)/*
#-----------------------------------------------------------------
# dependencies
#-----------------------------------------------------------------

depend: $(BLDSRC:.F90=.f90)
	~/bin/mkdepend-f90.awk $^ > $@.inc

-include depend.inc

depend_tl: $(TLMSRC:.F90=.f90)
	~/bin/mkdepend-f90.awk $^ > $@.inc

-include depend_tl.inc

depend_ad: $(ADMSRC:.F90=.f90)
	~/bin/mkdepend-f90.awk $^ > $@.inc

-include depend_ad.inc

depend_dot: $(DOTSRC:.F90=.f90)
	~/bin/mkdepend-f90.awk $^ > $@.inc

-include depend_dot.inc
