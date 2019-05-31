#!/bin/bash

export MKL_NUM_THREADS=1

for (( k=100; k<=1000; k+=100 ))
do
./xmain_A2QR_lapack.exe -m 8000 -n 8000 -k $k
done
for (( k=2000; k<=8000; k+=1000 ))
do
./xmain_A2QR_lapack.exe -m 8000 -n 8000 -k $k
done
