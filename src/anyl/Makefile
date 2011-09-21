#=======================================================================
#
# FILE:        Makefile
#
# SUMMARY:     Configurable Makefile for Enzo analysis toolkit 
#		(anyl, previously also known as enzo_anyl)
#
# AUTHOR:      Ji-hoon Kim (me@jihoonkim.org)
#
# DATE:        2009-07-29
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
EXE      = anyl
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

#MAKE_CONFIG_MACHINE  = $(ENZO_DIR)/Make.config.machine
include $(ENZO_DIR)/Make.config.machine

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing Enzo with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(ENZO_DIR)/Make.config.override
include $(MAKE_CONFIG_OVERRIDE)

# These may be necessary, so we manually override them here. --JHK
# No, we do not want to override defaults. --CBH
#CONFIG_USE_MPI   = yes
#CONFIG_PHOTON    = yes
#CONFIG_FAST_SIB  = yes
#CONFIG_LCAPERF   = no

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

$(EXE).exe: $(MODULES) dep $(OBJS_ANYL_LIB) purge_main
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) -o $(EXE).exe $(OBJS_ANYL_LIB) $(ENZO_DIR)/*.o \
             $(ENZO_DIR)/hydro_rk/*.o $(ENZO_DIR)/uuid/*.o $(LIBS) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
         else \
             echo "$(LD) $(LDFLAGS) -o $(EXE).exe $(OBJS_ANYL_LIB) $(ENZO_DIR)/*.o $(ENZO_DIR)/hydro_rk/*.o $(ENZO_DIR)/uuid/*.o $(LIBS)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

.PHONY: purge_main
purge_main:
	if [ -e $(ENZO_DIR)/enzo.o ]; then \
		/bin/rm $(ENZO_DIR)/enzo.o; \
		fi
	if [ -e $(ENZO_DIR)/EvolveHierarchy.o ]; then \
		/bin/rm $(ENZO_DIR)/EvolveHierarchy.o; \
		fi

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
             echo "$(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.c"; \
             echo; \
             $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.c;\
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
	-$(MAKE) show-flags >& $(TOP_DIR)/bin/$(EXE).show-flags
	-$(MAKE) show-config >& $(TOP_DIR)/bin/$(EXE).show-config
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
	@echo "   gmake                Compile and generate the executable 'anyl.exe'"
	@echo "   gmake install        Copy the executable to bin"
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
