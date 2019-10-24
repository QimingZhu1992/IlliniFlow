#!/bin/bash

FENICS_PYTHON=python3.5
TAG=2019.1.0
FENICS_ROOT_DIR=/home/zqm/fenics
PREFIX=${FENICS_ROOT_DIR}/fenics-${TAG}
BUILD_DIR=${FENICS_ROOT_DIR}/fenics-build

export PATH=${PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${PREFIX}/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${PREFIX}/include:${CPLUS_INCLUDE_PATH}

cd ${BUILD_DIR}
tar -xvf zlib-1.2.8.tar.gz
mkdir zlib
cd zlib-1.2.8
./configure --prefix=${BUILD_DIR}/zlib
make && make install
export PATH=${BUILD_DIR}/zlib/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/zlib/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/zlib/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/zlib/include:${CPLUS_INCLUDE_PATH}

cd ${BUILD_DIR}
tar -xvf boost_1_63_0.tar.gz
mkdir boost
cd boost_1_63_0
./bootstrap.sh --prefix=${BUILD_DIR}/boost
./b2
./b2 install
export PATH=${BUILD_DIR}/boost/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/boost/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/boost/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/boost/include:${CPLUS_INCLUDE_PATH}

cd ${BUILD_DIR}
tar -xvf pcre-8.42.tar.gz
mkdir pcre
cd pcre-8.42
./configure --prefix=${BUILD_DIR}/pcre
make
make install
export PATH=${BUILD_DIR}/pcre/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/pcre/lib:${LD_LIBRARY_PATH}
export LD_RUN_PATH=${BUILD_DIR}/pcre/lib:${LD_RUN_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/pcre/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/pcre/include:${CPLUS_INCLUDE_PATH}

cd ${BUILD_DIR}
tar -xvf pkg-config-0.29.2.tar.gz
tar -xvf bison-3.2.1.tar.gz
tar -xvf flex-2.6.4.tar.gz
mkdir pkg-config
mkdir bison
mkdir flex

cd ${BUILD_DIR}/pkg-config-0.29.2
./configure --prefix=${BUILD_DIR}/pkg-config
make && make install
export PATH=${BUILD_DIR}/pkg-config/bin:${PATH}
cd ..

cd ${BUILD_DIR}/bison-3.2.1
./configure --prefix=${BUILD_DIR}/bison
make && make install
export PATH=${BUILD_DIR}/bison/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/bison/lib:${LD_LIBRARY_PATH}
cd ..

cd ${BUILD_DIR}/flex-2.6.4
./configure --prefix=${BUILD_DIR}/flex
make && make install
export PATH=${BUILD_DIR}/flex/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/flex/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/flex/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/flex/include:${CPLUS_INCLUDE_PATH}
cd ..

cd ${BUILD_DIR}
mkdir cmake
tar -xvf cmake-3.11.0.tar.gz
cd ${BUILD_DIR}/cmake-3.11.0
./configure --prefix=${BUILD_DIR}/cmake
make && make install
export PATH=${BUILD_DIR}/cmake/bin:${PATH}
export LD_LIBRARY_PATH=${BUILD_DIR}/cmake/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${BUILD_DIR}/cmake/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${BUILD_DIR}/cmake/include:${CPLUS_INCLUDE_PATH}



