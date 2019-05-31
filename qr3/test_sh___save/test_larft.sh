#!/bin/bash

export MKL_NUM_THREADS=1

for (( m=1000; m<=30000; m+=1000 ))
do
./xmain_study_of_larft__lapack.exe -m $m -n 200 -k 200
done

for (( m=1000; m<=30000; m+=1000 ))
do
./xmain_study_of_larft__new.exe -m $m -n 200 -k 200
done

for (( m=1000; m<=30000; m+=1000 ))
do
./xmain_study_of_larft__qr3.exe -m $m -n 200 -k 200
done
