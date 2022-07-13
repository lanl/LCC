#!/bin/bash

# Make sure all the paths are correct

MY_PATH=$(pwd)
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export BLAS_VENDOR=${BLAS_VENDOR:=Generic}
export BML_OPENMP=${BML_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/bml/install"}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}

cd bml; ./build.sh configure; cd build; make -j; make install cd $MY_PATH


