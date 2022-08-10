#!/bin/bash

# Make sure all the paths are correct

MY_PATH=$(pwd)
BML_LIB="${MY_PATH}/bml/install"
export PKG_CONFIG_PATH="$BML_LIB/lib/pkgconfig:$BML_LIB/lib64/pkgconfig"
#export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:="$BML_LIB/lib/:$BML_LIB/lib64/"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export BLAS_VENDOR=${BLAS_VENDOR:=None}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/qmd-progress/install"}
export PROGRESS_TESTING=${PROGRESS_TESTING:=no}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=no}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}

cd qmd-progress; ./build.sh configure; cd build; make -j ;make install; cd $MY_PATH

                                                                                                                                                                                              
                                    
