#!/bin/bash

SALIDA=salida

export GFORTRAN_UNBUFFERED_ALL=1
export OMP_NUM_THREADS=6
export LIOHOME=/home/gonzalo/progs/LIOs/LR-fitting
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber

${LIOHOME}/liosolo/liosolo -i diazi.in -c diazi.xyz -v > $SALIDA
