#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

export GFORTRAN_UNBUFFERED_ALL=1
export OMP_NUM_THREADS=1
export LIOHOME=/home/gonzalo/progs/LIOs/openLR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOHOME/g2g:$LIOHOME/lioamber
${LIOHOME}/liosolo/liosolo -i li.in -c li.xyz -v > $SALIDA
