#!/bin/bash

FENICS_PYTHON=python3.5
TAG=2019.1.0
FENICS_ROOT_DIR=/home/zqm/fenics
PREFIX=${FENICS_ROOT_DIR}/fenics-${TAG}
BUILD_DIR=${FENICS_ROOT_DIR}/fenics-build
FENICS_VERSION="2019.1.0"

export PATH=${PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${PREFIX}/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${PREFIX}/include:${CPLUS_INCLUDE_PATH}

export PATH=${BUILD_DIR}/zlib/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/zlib/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/zlib/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/zlib/include:${CPLUS_INCLUDE_PATH}

export PATH=${BUILD_DIR}/boost/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/boost/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/boost/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/boost/include:${CPLUS_INCLUDE_PATH}
export BOOST_ROOT=${BUILD_DIR}/boost

export PATH=${BUILD_DIR}/pcre/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/pcre/lib:${LD_LIBRARY_PATH}
export LD_RUN_PATH=${BUILD_DIR}/pcre/lib:${LD_RUN_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/pcre/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/pcre/include:${CPLUS_INCLUDE_PATH}

export PATH=${BUILD_DIR}/pkg-config/bin:${PATH}

export PATH=${BUILD_DIR}/bison/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/bison/lib:${LD_LIBRARY_PATH}

export PATH=${BUILD_DIR}/flex/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/flex/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/flex/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/flex/include:${CPLUS_INCLUDE_PATH}

export PATH=${BUILD_DIR}/cmake/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/cmake/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/cmake/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/cmake/include:${CPLUS_INCLUDE_PATH}

export CC=mpicc
export CXX=mpic++

export PETSC_DIR=${PREFIX}

source ${PREFIX}/share/dolfin/dolfin.conf

