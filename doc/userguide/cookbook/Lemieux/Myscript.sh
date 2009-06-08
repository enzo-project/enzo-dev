#!/bin/csh
#PBS -l walltime=5:00:00
#PBS -l rmsnodes=4:16
#PBS -j oe
set echo
# execute program
prun -N ${RMS_NODES} -n ${RMS_PROCS} ./a.out

#  "walltime" tells the machine how much wall clock time is needed ofr
#     the job.
#  "rmsnodes"=(number of nodes):(total number of processors)
#  "-j oe" indicates that stderr and stdout will be writeen into a file.  That 
#     file is usually called something like MyScript.sh.o234, 
#     where 234 is the job number.

