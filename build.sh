# install conda and python dependencies
wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh -b -p ./enzo-conda -f
export PATH=$PWD/enzo-conda/bin:$PATH
conda install -q -y mercurial cython h5py matplotlib sympy numpy pytest flake8 yt nose
pip install python-hglib

# install OS dependencies
sudo apt-get update
sudo apt-get install -y csh libhdf5-serial-dev gfortran libtool openmpi-bin libopenmpi-dev

# install hypre
wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.11.2.tar.gz
tar xvfz hypre-2.11.2.tar.gz
cd hypre-2.11.2/src
./configure --prefix=/usr --with-MPI --with-MPI-include=/usr/include/mpi --with-MPI-libs=mpi --with-MPI-lib-dirs=/usr/lib
make install
cd ../../

# install grackle
mkdir -p $HOME/local
hg clone https://bitbucket.org/grackle/grackle
cd grackle
hg up tip
./configure
cd src/clib
make machine-linux-gnu
make
make install
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

echo "backend : Agg" > $HOME/matplotlibrc
export MATPLOTLIBRC=$HOME

cd $BITBUCKET_CLONE_DIR
hg up tip

./configure
cd src/enzo
make machine-linux-gnu
make opt-high
make integers-32
make particle-id-32
make hypre-yes
make grackle-yes
make -j 4
