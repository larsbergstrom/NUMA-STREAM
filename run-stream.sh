#!/bin/bash

for i in 1 2 4 8 12 24 36 48
do
  gcc -O2 -std=c99 -lpthread -lnuma -L ../numactl-2.0.7-rc2 -DNUM_THREADS=$i -DMEM_OFF=0 stream-pthreads.c
  echo "No off $i"
  ./a.out > NO_OFF_$i
done


for i in 1 2 4 8 12 24 36 48
do
  gcc -O2 -std=c99 -lpthread -lnuma -L ../numactl-2.0.7-rc2 -DNUM_THREADS=$i -DMEM_OFF=1 stream-pthreads.c
  echo "Off $i"
  ./a.out > OFF_$i
done

for i in 1 2 4 8 12 24 36 48
do
  gcc -O2 -std=c99 -lpthread -lnuma -L ../numactl-2.0.7-rc2 -DNUM_THREADS=$i -DMEM_OFF=0 -DSTRIDE=1 stream-pthreads.c
  echo "No off $i"
  ./a.out > NO_OFF_NO_STRIDE_$i
done


for i in 1 2 4 8 12 24 36 48
do
  gcc -O2 -std=c99 -lpthread -lnuma -L ../numactl-2.0.7-rc2 -DNUM_THREADS=$i -DMEM_OFF=1 -DSTRIDE=1 stream-pthreads.c
  echo "Off $i"
  ./a.out > OFF_NO_STRIDE_$i
done