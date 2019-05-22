#!/bin/bash

#for (( n=4; n<=200; n++ ))
#do
#./xmain_flops_geqr3.exe -m $n -n $n 
#done

./xmain_flops_geqr3_ISW.exe -m 300 -n 4
./xmain_flops_geqr3_ISW.exe -m 300 -n 8
./xmain_flops_geqr3_ISW.exe -m 300 -n 16
./xmain_flops_geqr3_ISW.exe -m 300 -n 32
./xmain_flops_geqr3_ISW.exe -m 300 -n 64
./xmain_flops_geqr3_ISW.exe -m 300 -n 128
./xmain_flops_geqr3_ISW.exe -m 300 -n 256

./xmain_flops_geqr3_ISW.exe -m 20 -n 4
./xmain_flops_geqr3_ISW.exe -m 20 -n 8
./xmain_flops_geqr3_ISW.exe -m 20 -n 16

./xmain_flops_geqr3_ISW.exe -m 40 -n 4
./xmain_flops_geqr3_ISW.exe -m 40 -n 8
./xmain_flops_geqr3_ISW.exe -m 40 -n 16
./xmain_flops_geqr3_ISW.exe -m 40 -n 32

./xmain_flops_geqr3_ISW.exe -m 70 -n 4
./xmain_flops_geqr3_ISW.exe -m 70 -n 8
./xmain_flops_geqr3_ISW.exe -m 70 -n 16
./xmain_flops_geqr3_ISW.exe -m 70 -n 32
./xmain_flops_geqr3_ISW.exe -m 70 -n 64

./xmain_flops_geqr3_ISW.exe -m 130 -n 4
./xmain_flops_geqr3_ISW.exe -m 130 -n 8
./xmain_flops_geqr3_ISW.exe -m 130 -n 16
./xmain_flops_geqr3_ISW.exe -m 130 -n 32
./xmain_flops_geqr3_ISW.exe -m 130 -n 64
./xmain_flops_geqr3_ISW.exe -m 130 -n 128





