#!/bin/sh
FFTW_DIR=/usr/work/jwise/local
./configure \
    CFLAGS="-I ${FFTW_DIR}/include" \
    FCFLAGS="-I ${FFTW_DIR}/include" \
    LDFLAGS="-L ${FFTW_DIR}/lib -lmpi" \
    FC=ifort \
    CC=icc \
    --enable-enzo \
    --enable-double \
    --with-hdf=/usr/work/jwise/local/hdf5/1.8.2p
