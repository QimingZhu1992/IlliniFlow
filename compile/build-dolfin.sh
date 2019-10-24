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

export CC=mpicc
export CXX=mpic++

cd ${BUILD_DIR}
VERSION="3.9.4"
tar -xf petsc-${VERSION}.tar.gz
cd petsc-${VERSION}
python2 ./configure --with-cc=mpicc \
        --with-fc=mpif90 \
        --with-cxx=mpicxx \
        --COPTFLAGS="-O2" \
        --CXXOPTFLAGS="-O2" \
        --FOPTFLAGS="-O2" \
	--download-cmake \
	--download-fblaslapack=1 \
        --download-hdf5 \
        --download-metis \
        --download-parmetis \
        --download-suitesparse \
        --download-scalapack \
        --download-hypre \
        --download-mumps \
        --download-ml \
        --download-scotch \
        --download-ptscotch \
        --with-debugging=0 \
        --with-shared-libraries \
        --prefix=${PREFIX}
make && make install
cd ..

export PETSC_DIR=${PREFIX}
export PETSC_ARCH=

cd ${BUILD_DIR}
SWIG_VERSION="3.0.12"
tar xf swig-${SWIG_VERSION}.tar.gz
cd swig-${SWIG_VERSION}
./configure --prefix=${PREFIX}
make && make install
cd ..

cd ${BUILD_DIR}
EIGEN_VERSION="3.3.3"
mkdir -p ${BUILD_DIR}/eigen
tar -xf eigen.tar.bz2 -C ${BUILD_DIR}/eigen --strip-components=1
cd eigen
mkdir -p build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=${PREFIX}
make install
cd ..

cd ${BUILD_DIR}
tar -xvf numpy-1.17.3.tar.gz
cd numpy-1.17.3
${FENICS_PYTHON} setup.py install --user

cd ${BUILD_DIR}
tar -xvf mpmath-1.1.0.tar.gz
cd mpmath-1.1.0
${FENICS_PYTHON} setup.py install --user

cd ${BUILD_DIR}
tar -xvf sympy-1.1.1.tar.gz
cd sympy-1.1.1
${FENICS_PYTHON} setup.py install --user

cd ${BUILD_DIR}
tar -xvf ply-3.11.tar.gz
cd ply-3.11
${FENICS_PYTHON} setup.py install --user

cd ${BUILD_DIR}
tar -xvf petsc4py-3.9.1.tar.gz
cd petsc4py-3.9.1
${FENICS_PYTHON} setup.py install --user
cd ..

cd $BUILD_DIR && \
    git clone https://bitbucket.org/fenics-project/dijitso.git && \
    cd dijitso && \
    git checkout ${FENICS_VERSION} && \
    ${FENICS_PYTHON} setup.py install --user

cd $BUILD_DIR && \
   git clone https://bitbucket.org/fenics-project/ufl.git
   cd ufl && \
   git checkout ${FENICS_VERSION} && \
   ${FENICS_PYTHON} setup.py install --user

cd $BUILD_DIR && \ 
   git clone https://bitbucket.org/fenics-project/fiat.git && \
   cd fiat && \
   git checkout ${FENICS_VERSION} && \
   ${FENICS_PYTHON} setup.py install --user

cd $BUILD_DIR && \
    git clone https://bitbucket.org/fenics-project/ffc.git && \
    cd ffc && \
    git checkout ${FENICS_VERSION} && \
    ${FENICS_PYTHON} setup.py install --user

USE_PYTHON3=on
cd $BUILD_DIR && \
    git clone https://bitbucket.org/fenics-project/dolfin.git && \
    cd dolfin && \
    git checkout ${FENICS_VERSION} && \
    mkdir -p build && \
    cd build && \
    cmake ../ -DDOLFIN_ENABLE_DOCS=False -DDOLFIN_ENABLE_OPENMP=off -DDOLFIN_ENABLE_MPI=on -DDOLFIN_SKIP_BUILD_TESTS=True -DCMAKE_BUILD_TYPE=Release -DSLEPC_INCLUDE_DIRS=${PREFIX}/include -DPETSC_INCLUDE_DIRS=${PREFIX}/include -DSWIG_EXECUTABLE:FILEPATH=${PREFIX}/bin/swig -DEIGEN3_INCLUDE_DIR:FILEPATH=${PREFIX}/include/eigen3 -DCMAKE_INSTALL_PREFIX=${PREFIX} -DPYTHON_EXECUTABLE:FILEPATH=$(which ${FENICS_PYTHON}) && \
    make && make install && \
    cd ..

source ${PREFIX}/share/dolfin/dolfin.conf
