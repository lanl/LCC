#!/bin/bash

set -e -u -x

SUDO=$(which sudo || true)

for i in $(seq 5); do
  ${SUDO} apt-get update && break
done

${SUDO} apt-get install --assume-yes --no-install-recommends \
  build-essential \
  bundler \
  cmake cmake-data \
  gcc-9 g++-9 gfortran-9 \
  gcc-10 g++-10 gfortran-10 \
  gcc-11 g++-11 gfortran-11 \
  git-core \
  indent \
  libblas-dev \
  liblapack-dev \
  libmetis-dev \
  libopenmpi-dev \
  make \
  pkg-config \
  python3-numpy \
  sudo
