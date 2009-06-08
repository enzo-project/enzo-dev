#!/bin/sh

#################### User-defined parameters ####################
PPN=2                     		   # Processes per node
MODULES="comp/intel-10.1.021 mpi/impi-3.2.011"  # Modules to load
OUTNAME="estd.out"  			   # Basename for STDOUT/STDERR redirect
PROCYON=1   				   # To use procyon for file management
WALLTIME="12:00:00"                        # Walltime to request
EXE=../enzo.exe
MACHINE=discover
#################################################################

if [ $# -lt 2 ]; then
    echo "usage: $0 ncpus amr_file"
    exit 1;
fi

NCPUS=$1
AMRFILE=$2
BASEDIR=`pwd`

# If the files have been moved to a directory, move them back
#if [ -d ${AMRFILE}.dir ]; then
#    if [ -f ${AMRFILE}.dir/${AMRFILE} ]; then
#	cd ${AMRFILE}.dir
#	echo "moving ${AMRFILE} data to ${BASEDIR}"
#	find . -maxdepth 1 -type -f -name "${AMRFILE}*" | \
#	    xargs mv --target-directory=..
#    fi
#fi

if [ ! -f ${AMRFILE} ]; then
    echo "${AMRFILE} not found.  Not re-submitting."
    exit 1;
fi

# Extract data dump number
NCHARS=`echo $AMRFILE | wc -c | gawk '{printf "%d", $1}'`
C1=`expr $NCHARS - 4`
C2=`expr $NCHARS - 5`
NUM=`echo $AMRFILE | cut -c${C1}- | gawk '{printf "%d", $1}'`
PREFIX=`echo $AMRFILE | cut -c1-${C2}`

if [ $NCPUS -lt $PPN ]; then
    PPN=$NCPUS
fi

QSUB=reqsub.sh
JOB_NAME=${BASEDIR##*/}
NODES=`expr $NCPUS / $PPN`
OUTNAME=${OUTNAME}.${NUM}

# discover at NCCS
if [ "$MACHINE" = "discover" ]; then
    if [ $NCPUS -le 16 ]; then
	QUEUE=general_small
    else
	QUEUE=general
    fi
    SELECT="#PBS -l select=${NODES}:ncpus=${PPN},walltime=${WALLTIME}"
fi

# palm at NCCS
if [ "$MACHINE" = "palm" ]; then
    if [ $NCPUS -le 16 ]; then
	QUEUE=general_small
    else
	QUEUE=general
    fi
    SELECT="#PBS -l ncpus=${NCPUS},walltime=${WALLTIME}"
fi

if [ $PROCYON -gt 0 ]; then
    PSCRIPT="cp -f ../procyon . && ./procyon -600 &"
else
    PSCRIPT=""
fi

rm -f $QSUB
cat >> $QSUB << EOF
#!/usr/local/bin/csh
#PBS -N ${JOB_NAME}
$SELECT
#PBS -l place=scatter
#PBS -j oe
#PBS -m be
#PBS -q $QUEUE
#PBS -r n
#PBS -W group_list=g0760
#PBS -v SCAFUN_CACHING_MODE=0,MPI_REQUEST_MAX=65536
module add $MODULES
cd \$PBS_O_WORKDIR
$PSCRIPT
mpirun -np $NCPUS $EXE -d -r $AMRFILE >& $OUTNAME
exit 0
EOF

qsub $QSUB