#!/bin/sh
export OPENBLAS_NUM_THREADS=1
for ((i=100; i<1000 ; i+=100));
do
./xmain_wLAPACK.exe -m $i -n $i -ii 0
done
#for ((i=2000; i<6000 ; i+=1000));
#do
#./xmain_wLAPACK.exe -m $i -n $i -ii 0
#done
