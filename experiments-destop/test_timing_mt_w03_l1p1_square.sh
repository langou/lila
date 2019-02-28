#!/bin/sh
export OPENBLAS_NUM_THREADS=1
rm -f test_timing_w03_l1p1__square__1thread.dat
rm -f test_timing_w03_l1p1__square__2threads.dat
rm -f test_timing_w03_l1p1__square__3threads.dat
rm -f test_timing_w03_l1p1__square__4threads.dat
for ((i=1; i<1100 ; i+=100));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__1thread.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__1thread.dat
done

export OPENBLAS_NUM_THREADS=2
for ((i=1; i<1100 ; i+=100));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__2threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__2threads.dat
done

export OPENBLAS_NUM_THREADS=3
for ((i=1; i<1100 ; i+=100));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__3threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__3threads.dat
done

export OPENBLAS_NUM_THREADS=4
for ((i=1; i<1100 ; i+=100));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__4threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_w03.exe -m 5000 -n 5000 -mt $i -nx 20 -ii 0 -leaf 1 -panel 1 >> test_timing_w03_l1p1__square__4threads.dat
done
