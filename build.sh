# install conda and python dependencies
wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh -b -p ./enzo-conda -f
export PATH=$PWD/enzo-conda/bin:$PATH
conda install -q -y mercurial cython h5py matplotlib sympy numpy pytest flake8 yt

# install OS dependencies
sudo apt-get update
sudo apt-get install -y csh libhdf5-serial-dev gfortran libtool openmpi-bin libopenmpi-dev

echo "backend : Agg" > $HOME/matplotlibrc
export MATPLOTLIBRC=$HOME

cd $BITBUCKET_CLONE_DIR
hg up tip

./configure
cd src/enzo
make machine-linux-gnu
make

