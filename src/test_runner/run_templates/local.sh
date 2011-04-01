#! /bin/bash

mpirun ./${ENZO_EXE} -n ${N_PROCS} ${PAR_FILE} >& estd.out;