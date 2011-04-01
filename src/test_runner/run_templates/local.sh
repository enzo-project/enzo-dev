#! /bin/bash

mpirun -n ${N_PROCS} ${ENZO_EXE} -d ${PAR_FILE} >& estd.out;
