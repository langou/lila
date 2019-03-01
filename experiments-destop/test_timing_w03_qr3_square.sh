#!/bin/sh

export OPENBLAS_NUM_THREADS=1

rm -f test_timing_qr3__square__1thread.dat
rm -f test_timing_w03_qr3__square__1thread.dat

for ((i=100; i<1100 ; i+=100));
do
../lila-c/xmain_w03_qr3.exe -m $i -n $i -mt 200 -ii 0 >> test_timing_w03_qr3__square__1thread.dat
done

for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_w03_qr3.exe -m $i -n $i -mt 200 -ii 0 >> test_timing_w03_qr3__square__1thread.dat
done


