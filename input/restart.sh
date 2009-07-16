#!/bin/sh

NPROCS=$1
AMRFILE=$2
EXE=/Users/jwise/codes/EnzoWOC/week-of-code/src/enzo/enzo.exe

# Extract data dump number
NCHARS=`echo $AMRFILE | wc -c | gawk '{printf "%d", $1}'`
C1=`expr $NCHARS - 4`
C2=`expr $NCHARS - 5`
NUM=`echo $AMRFILE | cut -c${C1}- | gawk '{printf "%d", $1}'`
PREFIX=`echo $AMRFILE | cut -c1-${C2}`

pwd
date

echo "mpirun -np $NPROCS $EXE -d -r $AMRFILE"
mpirun -np $NPROCS $EXE -d -r $AMRFILE &
