#!/bin/sh

rm -f test.dat

# export OPENBLAS_NUM_THREADS=18
export MKL_NUM_THREADS=1

for ((i=100; i<1100 ; i+=100));
do
for ((j=0; j<40 ; j++));
do
./xmain_dgeqrf.exe -m $i -n $i -k $i >> test.dat
done
done

for ((i=2000; i<10000 ; i+=500));
do
for ((j=0; j<10 ; j++));
do
./xmain_dgeqrf.exe -m $i -n $i -k $i >> test.dat
done
done

