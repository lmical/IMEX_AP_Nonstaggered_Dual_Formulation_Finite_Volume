#==================================================================================================#
#  MAKEFILE
#==================================================================================================#
#
#--------------------------------------------------------------------------------------------------#
INFO_COMPILING=\e[1;34m[COMPILING]\e[0m
INFO_LINKING=\e[1;34m[LINKING]\e[0m
#-------------------------------------------------------------------------------------------------#
COMPILER=GNU DEBUG
#--------------------------------------------------------------------------------------------------#
ifeq ($(findstring GNU,$(COMPILER)),GNU)
  FC=gfortran
  INSTRUCTION=-march=native
  FCFLAGS  = -xf95-cpp-input -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-300 -ffpe-trap=invalid  -ffree-line-length-none
  FLFLAGS  = -fdefault-real-8 -fbackslash
  ifeq ($(findstring DEBUG,$(COMPILER)),DEBUG)
    FCFLAGS  += -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace
    FLFLAGS  += -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace
  else
    FCFLAGS  += -O3 $(INSTRUCTION) -pg
    FLFLAGS  += -O3 $(INSTRUCTION) -pg
  endif
endif
#--------------------------------------------------------------------------------------------------#
ifeq ($(findstring INTEL,$(COMPILER)),INTEL)
  FC=ifort
  INSTRUCTION=-O3 -xCORE-AVX2 -ip
  FCFLAGS  =-fpp -allow nofpp_comments -assume bscc -r8 -i4 -traceback -warn all
  FLFLAGS  =-r8 -i4 -traceback -assume bscc
  ifeq ($(findstring DEBUG,$(COMPILER)),DEBUG)
    FCFLAGS  +=-g -O0 -fpe0 -traceback
    FCFLAGS  +=-check all,noarg_temp_created,noformat,nooutput_conversion,pointer,bounds,uninit
    FCFLAGS  +=-init=snan -init=arrays
    FLFLAGS  +=-g -O0
  else
    FCFLAGS  +=$(INSTRUCTION)
    FLFLAGS  +=$(INSTRUCTION)
  endif
endif
#--------------------------------------------------------------------------------------------------#
AR       =ar
AR_FLAGS =cr
#--------------------------------------------------------------------------------------------------#
BLAS_LAPACK = -llapack -lblas
#--------------------------------------------------------------------------------------------------#
MAIN_LIB = libmain.a
#--------------------------------------------------------------------------------------------------#
# FCFLAGS += -I. -DTRICKSTOAVOIDDIVISIONSBYZERO
# FCFLAGS += -I. -DACTIVEFLUX
FCFLAGS += -I. -DCENTEREDPRIMITIVE
# FCFLAGS += -I. -DRECONSTRUCTFROMCONSERVED
# FCFLAGS += -I. -DDONOTOVERWRITEPRIMITIVE
# DO NOT USE IT # FCFLAGS += -I. -DEVOLVERHOWITHSAMESCHEME
# DO NOT USE IT # FCFLAGS += -I. -DOLDEVOLVERHOWITHSAMESCHEME
# FCFLAGS += -I. -DPATHCONSERVATIVESHOCKDETECTION
# FCFLAGS += -I. -DSMOOTHINGLIMITER
# FCFLAGS += -I. -DDETECTEACHSTAGE
# FCFLAGS += -I. -DFIRSTORDERPRIMITIVE
FCFLAGS += -I. -DRELAXATION
FCFLAGS += -I. -DIMEX
# FCFLAGS += -I. -DBLENDINGDT
# FCFLAGS += -I. -DNOMODIFICATIONWAVESPEEDIMEX
FCFLAGS += -I. -DHYPERBOLICITYTRICKMODIFIED
# FCFLAGS += -I. -DCENTEREDPRIMITIVEIMEX
# FCFLAGS += -I. -DIMEXMOMENTUM
# FCFLAGS += -I. -DPRIMITIVEONLY
# FCFLAGS += -I. -DMOMENTUMINPRIMITIVEVARIABLES
# FCFLAGS += -I. -DNORMALIZELINEARSYSTEMS
FCFLAGS += -I. -DNORMALIZEALL
# FCFLAGS += -I. -DMULTIFLUID
# FCFLAGS += -I. -DPWCINTROUBLEDCELLS
# FCFLAGS += -I. -DALTERNATIVEFORMULATIONPRESSURE
# FCFLAGS += -I. -DPRESSUREFORMULATIONFORIMEXMOMENTUM
# FCFLAGS += -I. -DRECONSTRUCTIONPHI
# FCFLAGS += -I. -DRECONSTRUCTIONVELOCITY
# FCFLAGS += -I. -DLEVELSETINU
# FCFLAGS += -I. -DPRIMITIVEFORMULATIONLEVELSET
# FCFLAGS += -I. -DOTHERPOSTPROCESSING
# FCFLAGS += -I. -DSIDEDPOSTPROCESSING
# FCFLAGS += -I. -DALINAPOSTPROCESSING
# FCFLAGS += -I. -DFOURINPUTSPOSTPROCESSING
# FCFLAGS += -I. -DOVERCOMPRESSIVEPOSTPROCESSING
# FCFLAGS += -I. -DPWLFORGAMMAPINF
# FCFLAGS += -I. -DNEWANTIDIFFUSIONTERM
# FCFLAGS += -I. -DPOSITIVITYCHECKS
# FCFLAGS += -I. -DPOSTPROCESSINGATEACHSTEP
# FCFLAGS += -I. -DCOUNTJACOBI
# FCFLAGS += -I. -DIMEXL2FULLYUPWINDED
# FCFLAGS += -I. -DSW
#FCFLAGS += -I. -DWELLBALANCED
#FCFLAGS += -I. -DPATANKAR
#--------------------------------------------------------------------------------------------------#
SRC = $(shell cat sources.mk)
OBJ = $(SRC:.f90=.o)
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
all:
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "COMPILING MAIN CODE...\n"
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@$(MAKE) -s $(MAIN_LIB)
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "$(INFO_COMPILING)   $(FC) -c main.f90 -o main.o\n";
	@$(FC) $(FCFLAGS) -c main.f90 -o main.o
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "LINKING MAIN CODE...\n"
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "$(INFO_LINKING)     $(FC) main.o ../lib/$(MAIN_LIB) -o ../bin/main\n";
	@$(FC) main.o  ../lib/$(MAIN_LIB) $(FLFLAGS) $(BLAS_LAPACK) -o ../bin/main
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "COMPILING AND LINKING MAIN CODE...DONE!\n"
	@printf "%`tput cols`s\n"|tr ' ' '-'
#--------------------------------------------------------------------------------------------------#
$(MAIN_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) ../lib/$(MAIN_LIB) $(OBJ) #$(FLFLAGS)
#--------------------------------------------------------------------------------------------------#
%.o: %.f90
	@printf "$(INFO_COMPILING)   $(FC) $<\n";
	@$(FC) $(FCFLAGS) -c $< -o $@
#--------------------------------------------------------------------------------------------------#
clean:
	rm -f $(OBJ) *.mod *.i main.o ../bin/main ../lib/$(MAIN_LIB)
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
# DEPENDENCIES
#--------------------------------------------------------------------------------------------------#
include deplist.mk
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
# Utility Targets
#--------------------------------------------------------------------------------------------------#
.PHONY: all clean 
#--------------------------------------------------------------------------------------------------#
