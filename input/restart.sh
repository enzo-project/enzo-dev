#!/bin/sh

# Initial parameter file.  Will change after restarts
AMRFILE=nBody.RunParamFile

# Set to 1 if the initial run is a restart
INIT_RESTART=0

# Redirection of output (doesn't work ...)
OUTPUT=""
#OUTPUT=">& estd.out"

# MPI call
MPI=mpirun

# Number of processors
NPROCS=1

# Enzo executable
EXE=/Users/jwise/codes/EnzoWOC/week-of-code/src/enzo/enzo.exe

# These shouldn't be changed as they're (for now) hard-coded into enzo
RESTART_FILE=RestartParamFile
FINISH_FILE=RunFinished

########################################################################
########################################################################

# Remove restarting and stopping flag files
if [ -e $RESTART_FILE ]; then
    rm -f $RESTART_FILE
fi
if [ -e $FINISH_FILE ]; then
    rm -f $FINISH_FILE
fi

if [ $INIT_RESTART -eq 1 ]; then
    RESTART_FLAG="-r"
else
    RESTART_FLAG=""
fi

# Loop until the run is finished.  Restarting when found the restart flag.
FIRST_TIME=1
while [ ! -f $FINISH_FILE ]; do
    echo "$MPI -np $NPROCS $EXE -d $RESTART_FLAG $AMRFILE $OUTPUT"
    $MPI -np $NPROCS $EXE -d $RESTART_FLAG $AMRFILE $OUTPUT

    # Exit if an error code is returned
    if [ "$?" -ne "0" ]; then
	exit 1;
    fi

    # Get restart parameter file from RESTART_FILE
    if [ -f $RESTART_FILE ]; then

	# If this isn't the first loop, erase the previous restart dump
	DIR=${AMRFILE%%/*}
	if [ -d $DIR -a $FIRST_TIME -eq 0 ]; then
	    rm -r $DIR
	fi

	AMRFILE=`cat $RESTART_FILE`
	RESTART_FLAG="-r"
	FIRST_TIME=0

    fi

done