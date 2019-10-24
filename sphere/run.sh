rm -rf CMakeFiles
rm -rf CMakeCache.txt
rm -rf cmake_install.cmake
rm -rf Makefile
rm -rf ns.h
rm -rf calcd.h
rm -rf calQ.h
ffc -l dolfin ns.ufl
ffc -l dolfin calcd.ufl
ffc -l dolfin calQ.ufl
cmake .
make

