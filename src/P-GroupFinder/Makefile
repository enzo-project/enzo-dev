#=======================================================================
#
# FILE:        Makefile
#
# SUMMARY:     Configurable Makefile for Enzo analysis toolkit 
#		(anyl, previously also known as enzo_anyl)
#
# AUTHOR:      John Wise (john.h.wise@nasa.gov)
#
# DATE:        2009-08-07
#
# DESCRIPTION
#	       Adopted from Makefile for enzohop
#	
#              See 'gmake help' for definitive description of targets
#
#              Makefile includes the following files:
# 
#              Make.config.settings   default configuration settings
#              Make.config.override   local user settings
#              Make.config.assemble   maps 'config' settings to 'flag' settings
#              Make.config.targets    configure targets
#              Make.mach.*            all machine-dependent settings
#              Make.config.objects    list of object files
#              DEPEND                 Make-generated dependencies
#
#              Make.mach.* should be the only file that one should
#              edit when porting Enzo to a new machine.
#
#              Make.config.override should be the only file that
#              one should edit when defining local user settings.
#              Preferably, this is done implicitly through
#              the available make targets (e.g. "gmake precision-32").
#              Use "gmake help-config" for a list of all configuration
#              settings.  These make targets do error-checking; hand-editing 
#              Make.config.override is more error-prone.
#
#=======================================================================

# Use bash since sh on datastar does not recognize ">&" used in dep: target

SHELL    = /bin/bash

TOP_DIR  = ../..
EXE      = P-GroupFinder
OUTPUT   = out.compile

ENZO_DIR = $(TOP_DIR)/src/enzo
MODULES  = 

#-----------------------------------------------------------------------
# Make.config.settings is used for setting default values to all compile-time 
# configuration settings.
#-----------------------------------------------------------------------

include $(ENZO_DIR)/Make.config.settings

#-----------------------------------------------------------------------
# Make.config.machine is used for setting which Make.mach.* file to use
#-----------------------------------------------------------------------

MAKE_CONFIG_MACHINE  = $(ENZO_DIR)/Make.config.machine
include $(ENZO_DIR)/Make.config.machine

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing Enzo with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(ENZO_DIR)/Make.config.override
include $(MAKE_CONFIG_OVERRIDE)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ANYL ONLY: override the override, These I found are needed 
# for smooth compiling, but you may want to try other options.

#CONFIG_PRECISION = 32
CONFIG_INTEGERS  = 32
CONFIG_USE_MPI   = yes
CONFIG_LCAPERF   = no
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include $(ENZO_DIR)/Make.config.assemble

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#-----------------------------------------------------------------------

-include $(ENZO_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.enzo/Make.mach.$(CONFIG_MACHINE)

#=======================================================================
# OBJECT FILES
#=======================================================================

include Make.config.objects

#-----------------------------------------------------------------------
# MAIN TARGET
#-----------------------------------------------------------------------

$(EXE): $(MODULES) dep $(OBJS_PGROUP_LIB)
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) -o $(EXE) $(OBJS_PGROUP_LIB) $(LIBS) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
         else \
             echo "$(LD) $(LDFLAGS) -o $(EXE) $(OBJS_PGROUP_LIB) $(LIBS)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

#-----------------------------------------------------------------------
# Ripped out the particle ID setup from targets/assemble because we
# ignore all of the other defines.
# -----------------------------------------------------------------------

VALID_IDS = particle-id-32 particle-id-64
.PHONY: $(VALID_IDS)

particle-id-32: CONFIG_PARTICLE_IDS-32
particle-id-64: CONFIG_PARTICLE_IDS-64
particle-id-%:
	@printf "\n\tInvalid target: $@\n\n\tValid targets: [$(VALID_IDS)]\n\n"
CONFIG_PARTICLE_IDS-%: suggest-clean
	@tmp=.config.temp; \
	grep -v CONFIG_PARTICLE_IDS $(MAKE_CONFIG_OVERRIDE) > $${tmp}; \
	mv $${tmp} $(MAKE_CONFIG_OVERRIDE); \
	echo "CONFIG_PARTICLE_IDS = $*" >> $(MAKE_CONFIG_OVERRIDE); \
	$(MAKE)  show-config | grep CONFIG_PARTICLE_IDS:
	@echo

#-----------------------------------------------------------------------
# DETERMINE INTEGER SIZE SETTINGS FOR PARTICLE INDICES
#-----------------------------------------------------------------------

ERROR_IDS = 1

# Settings for 32-bit integers

ifeq ($(CONFIG_PARTICLE_IDS),32)
    ERROR_IDS = 0
    ASSEMBLE_IDS_DEFINES = -DCONFIG_PINT_4
endif

# Settings for 64-bit integers

ifeq ($(CONFIG_PARTICLE_IDS),64)
    ERROR_IDS = 0
    ASSEMBLE_IDS_DEFINES = -DCONFIG_PINT_8
endif

# error if CONFIG_INTEGERS is incorrect

ifeq ($(ERROR_IDS),1)
   .PHONY: error_ids
   error_ids:
	$(error Illegal value '$(CONFIG_PARTICLE_IDS)' for $$(CONFIG_PARTICLE_IDS))
endif

#=======================================================================

# Override previous defines with the exception of particle id size
DEFINES = -DUSE_HDF5 -DFOF_ONLY -DH5_USE_16_API $(ASSEMBLE_IDS_DEFINES)

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .src .src90 .f .f90

.src.f:
	@echo "Compiling $<"
	@$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f
.f.o:
	@rm -f $@
	-@($(FC) -c -o $@ $(FFLAGS) $*.f) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             rm -f $*.f; \
          else \
             echo; \
             echo "$(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f"; \
             echo "$(FC) -c -o $@ $(FFLAGS) $*.f"; \
             echo; \
             $(CPP) $(DEFINES) $(CPPFLAGS) $*.src > $*.f; \
             $(FC) -c -o $@ $(FFLAGS) $*.f; \
             echo; \
             exit 1; \
          fi)

.src90.f90:
	@echo "Compiling $<"
	@$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90
.f90.o:
	@rm -f $@
	-@($(F90) -c -o $@ $(F90FLAGS) $*.f90) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             rm -f $*.f90; \
          else \
             echo; \
             echo "$(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90"; \
             echo "$(F90) -c -o $@ $(F90FLAGS) $*.f90"; \
             echo; \
             $(CPP) $(DEFINES) $(CPPFLAGS) $*.src90 > $*.f90; \
             $(F90) -c -o $@ $(F90FLAGS) $*.f90; \
             echo; \
             exit 1; \
          fi)

.c.o:
	@rm -f $@
	@echo "Compiling $<"
	-@($(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c) >& $(OUTPUT)
	@(if [ ! -e $@ ]; then \
             echo; \
             echo "$(CC) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.c"; \
             echo; \
             $(CC) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.c;\
             echo; \
             exit 1; \
          fi)

.C.o:
	@rm -f $@
	@echo "Compiling $<"
	-@($(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C) >& $(OUTPUT)
	@(if [ ! -e $@ ]; then \
             echo; \
             echo "$(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C"; \
             echo; \
             $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C;\
             echo; \
             exit 1; \
          fi)

#-----------------------------------------------------------------------
# Generate dependency file
#-----------------------------------------------------------------------

.PHONY: dep
dep:
	@echo "Updating DEPEND"
	-@(makedepend $(DEFINES) $(INCLUDES) -fDEPEND -m -- -- *.c *.h) >& out.make.DEPEND

include DEPEND

#-----------------------------------------------------------------------
# INSTALL TARGET
#-----------------------------------------------------------------------

install: $(TOP_DIR)/bin
	if [ ! -e $(TOP_DIR)/bin ]; then mkdir $(TOP_DIR)/bin; fi
	-@cp $(EXE) $(TOP_DIR)/bin/$(EXE)

$(TOP_DIR)/bin:
	mkdir $(TOP_DIR)/bin

#-----------------------------------------------------------------------
# HELP TARGET
#-----------------------------------------------------------------------

help:
	@echo
	@echo "========================================================================"
	@echo "   Enzo Makefile Help"
	@echo "========================================================================"
	@echo
	@echo "   gmake                Compile and generate the executable 'anyl.exe'"
	@echo "   gmake install        Copy the executable to bin"
	@echo "   gmake help           Display this help information"
	@echo "   gmake clean          Remove object files, executable, etc."
	@echo "   gmake dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   gmake help-config    Display detailed help on configuration make targets"
	@echo "   gmake show-config    Display the configuration settings"
	@echo "   gmake show-flags     Display specific compilation flags"
	@echo "   gmake default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------

clean:
	-@rm -f *.o *.mod *.f *.f90 nrsrc/*.o *~ $(OUTPUT) $(EXE)

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(ENZO_DIR)/Make.config.targets
