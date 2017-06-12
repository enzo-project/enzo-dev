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

export ENZOTEST_DIR=$HOME/enzo_test

# Build the gold standard version.
cd $BITBUCKET_CLONE_DIR
hg up test-gold-standard
./configure
cd src/enzo
make machine-linux-gnu
make load-config-test-suite
make -j 4

# Generate the gold standard results.
cd $BITBUCKET_CLONE_DIR/run
python ./test_runner.py --suite=push -o $ENZOTEST_DIR --answer-store --answer-name=push_suite  --local --strict=high --verbose

# Build the tip version.
cd $BITBUCKET_CLONE_DIR
hg up tip
./configure
cd src/enzo
make clean
make machine-linux-gnu
make load-config-test-suite
make -j 4

# Run tests on the tip and compare to gold standard.
cd $BITBUCKET_CLONE_DIR/run
python ./test_runner.py --suite=push -o $ENZOTEST_DIR --answer-name=push_suite  --local --clobber --strict=high --verbose
