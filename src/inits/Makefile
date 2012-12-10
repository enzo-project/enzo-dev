#=======================================================================
#
# FILE:        Makefile
#
# SUMMARY:     Configurable Makefile for Enzo inits
#
# AUTHOR:      James Bordner (jobordner@ucsd.edu)
#
# DATE:        2007-02-21
#
# DESCRIPTION
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
EXE      = inits
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
# INITS ONLY: override the override, since we don't want to compile with MPI

CONFIG_USE_MPI = no
CONFIG_LCAPERF = no

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

$(EXE).exe: $(MODULES) dep $(OBJS_LIB)
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) -o $(EXE).exe $(OBJS_LIB) $(LIBS) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
             if [ ! -e $(TOP_DIR)/bin ]; then mkdir $(TOP_DIR)/bin; fi; \
             cp $(EXE).exe $(TOP_DIR)/bin/$(EXE); \
         else \
             echo "$(LD) $(LDFLAGS) -o $(EXE).exe $(OBJS_LIB) $(LIBS)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .src .src90 .f .f90 .F .F90

.F.o:
	@echo "Compiling $<"
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(FC) -c -o $@ $(FFLAGS) $(DEFINES) $*.F >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(FC) -c -o $@ $(FFLAGS) $(DEFINES) $*.F"; \
             echo; \
             $(FC) -c -o $@ $(FFLAGS) $(DEFINES) $*.F; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(FC) -c -o $@ $(FFLAGS)  $(DEFINES) $*.F >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.F90.o:
	@rm -f $@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(F90) -c -o $@ $(F90FLAGS) $(DEFINES) $*.F90 >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(F90) -c -o $@ $(DEFINES) $(F90FLAGS) $*.F90"; \
             echo; \
             $(F90) -c -o $@ $(F90FLAGS) $*.F90; \
             echo; \
             exit 1; \
	  fi ; \
	else \
	  $(F90) -c -o $@ $(F90FLAGS) $(DEFINES) $*.F90 >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

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
             echo "$(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C"; \
             echo; \
             $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C;\
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
	-@(makedepend $(DEFINES) $(INCLUDES) -fDEPEND -m -- -- *.C *.src *.h) >& out.make.DEPEND

include DEPEND

#-----------------------------------------------------------------------
# INSTALL TARGET
#-----------------------------------------------------------------------

install: $(TOP_DIR)/bin
	if [ ! -e $(TOP_DIR)/bin ]; then mkdir $(TOP_DIR)/bin; fi
	-@cp $(EXE).exe $(TOP_DIR)/bin/$(EXE)
	-$(MAKE) show-flags   >& $(TOP_DIR)/bin/$(EXE).show-flags
	-$(MAKE) show-config  >& $(TOP_DIR)/bin/$(EXE).show-config
	-$(MAKE) show-version >& $(TOP_DIR)/bin/$(EXE).show-version

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
	@echo "   gmake                Compile and generate the executable 'enzo.exe'"
	@echo "   gmake install        Copy the executable to bin/enzo"
	@echo "   gmake help           Display this help information"
	@echo "   gmake clean          Remove object files, executable, etc."
	@echo "   gmake dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   gmake get-version    Generate version-related files"
	@echo "   gmake show-version   Display Subversion branch and revision"
	@echo
	@echo "   gmake help-config    Display detailed help on configuration make targets"
	@echo "   gmake show-config    Display the configuration settings"
	@echo "   gmake show-flags     Display specific compilation flags"
	@echo "   gmake default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------

clean:
	-@rm -f *.o *.mod *.f *.f90 *~ $(OUTPUT) $(EXE).exe

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(ENZO_DIR)/Make.config.targets
