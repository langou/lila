#!/bin/sh
export OPENBLAS_NUM_THREADS=1
rm -f test_timing_lapack__square__1thread.dat
rm -f test_timing_lapack__square__2threads.dat
rm -f test_timing_lapack__square__3threads.dat
rm -f test_timing_lapack__square__4threads.dat
for ((i=100; i<1100 ; i+=100));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__1thread.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__1thread.dat
done

export OPENBLAS_NUM_THREADS=2
for ((i=100; i<1000 ; i+=100));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__2threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__2threads.dat
done

export OPENBLAS_NUM_THREADS=3
for ((i=100; i<1000 ; i+=100));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__3threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__3threads.dat
done

export OPENBLAS_NUM_THREADS=4
for ((i=100; i<1000 ; i+=100));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__4threads.dat
done
for ((i=2000; i<6000 ; i+=1000));
do
../lila-c/xmain_wLAPACK.exe -m $i -n $i -ii 0 >> test_timing_lapack__square__4threads.dat
done
