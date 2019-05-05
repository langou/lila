#!/bin/bash

for (( n=10; n<=15; n++ ))
do
for (( nb=1; nb<=$n; nb++ ))
do
./xmain_flops_lapack_geqrf.exe -m 100 -n $n -nb $nb 
done
done
