#!/bin/sh
export OPENBLAS_NUM_THREADS=4
rm -f test_timing_lapack__square__4threads.dat
for ((i=100; i<1000 ; i+=100));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__4threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__4threads.dat
done
