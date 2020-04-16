#! /bin/bash
#
mpicxx -c -Wall hello_mpi.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mpicxx hello_mpi.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hello_mpi.o
mv a.out $HOME/bincpp/hello_mpi
#
echo "Normal end of execution."
