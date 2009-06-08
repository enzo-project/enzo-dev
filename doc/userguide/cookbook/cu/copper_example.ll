# Sample batch script for a mpi job on NCSA Copper
#
# Submit this script using the "llsubmit" command: 
#   llsubmit copper_example.ll
#
# Use the "llq" command to check the status of a job.
#
# Lines starting with #@ are embedded loadleveler directives.
# To comment out a loadleveler directive, put # in front of #@ 

#@ shell = /usr/bin/tcsh
#@ job_type = parallel

# copy all environment variables to batch job
#@ environment = COPY_ALL

# Send email when job completes
# (notification options: always,start,error,complete,never)
#@ notification = complete

# Specify job class
#@ class = batch

# Specify number of MPI processes
#@ tasks_per_node = 16

# Specify memory = amount per MPI process
# For mpi, ConsumableCpus is always 1, and it must be specified.
#@ resources = ConsumableCpus(1) ConsumableMemory(1000Mb)

# Specify the wall clock limit = hrs:min:sec
#@ wall_clock_limit = 12:00:00

# Specify the name of the job
#@ job_name = mpps_1

# Specify the standard output and standard error for the job
# They will be written in the subdirectory from which llsubmit is run
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err

# This has to be the last loadleveler directive
#@ queue
##########################################################


# change to the scratch directory
cd $SCR

# get executable and input file from user home directory

cp /u/ac/bwoshea/mpps/ics/* .
cp /u/ac/bwoshea/mpps/parfiles/mncp_parfile_starform.1 .

cp /u/ac/bwoshea/E8.6/amr_mpi/exe/enzo .
cp /u/ac/bwoshea/E8.6/amr_mpi/exe/cool_rates.in .
cp /u/ac/bwoshea/E8.6/amr_mpi/exe/ATOMIC.DAT .

# due to unitree "feature", must set executable bit
chmod 700 enzo

# back up to unitree (mss) if it crashes!
msscmd -b -f - << EOF

cd /u/ac/bwoshea/mpps_feb03/starform/1/
tar cvf sf_mss_1.tar *log *out perf*

EOF

#
# backup script - backs up data to unitree while
# enzo is running (the & forces it to the background)
squirrel /u/ac/bwoshea/mpps_feb03/starform/1/ &

# set limits on memory use and core size
limit coredumpsize 40000
limit memoryuse unlimited

# run mpi executable 
poe enzo -d mncp_parfile_starform.1 > mpps.1.log

# after run is over, put stuff into mass storage
tar cvf sf_1.tar *log *out perf*
msscmd cd /u/ac/bwoshea/mpps_feb03/starform/1/, put sf_1.tar

# cause backup script to crash
sleep 3600
touch squirrel_stop
sleep 100
