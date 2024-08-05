#PBS -S /bin/sh
#PBS -N enzo_compile
#PBS -l select=1:ncpus=28:model=bro
#PBS -l walltime=2:00:00
#PBS -q devel
#PBS -j oe
#PBS -o /home5/clochhaa/FOGGIE/output_compile_enzo
#PBS -koed
#PBS -m abe
#PBS -V
#PBS -W group_list=s2961
#PBS -l site=needed=/home5+/nobackupp13

# Be sure to change the above path to the output file and this cd to your own enzo-foggie directory!
cd /home5/clochhaa/enzo-foggie/src/enzo

module load comp-intel/2020.4.304
module load hdf5/1.8.18_serial

make machine-pleiades-mpich
make grackle-yes
make opt-high
make clean
make -j16
