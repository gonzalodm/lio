#!/bin/bash

SALIDA=salida

export GFORTRAN_UNBUFFERED_ALL=1
export OMP_NUM_THREADS=1
export LIOHOME=/home/gonzalo/progs/LIOs/LR-fitting
#export LIOHOME=/home/gonzalo/progs/LIOs/PCG_SOLVE
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber

#for num in {1..50}
#do
#  echo $num
  ${LIOHOME}/liosolo/liosolo -i hemo.in -c hemo.xyz -v > $SALIDA
#  grep "CONVERGED" salida
#  mv salida salida_$num
#done
