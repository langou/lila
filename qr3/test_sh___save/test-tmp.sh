#!/bin/bash
for (( m=11; m<=20; m++ ))
do
for (( k=5; k<=$m; k+=2 ))
do
for (( n=$k; n<=$m; n++ ))
do
for (( nb=1; nb<=10; nb+=3 ))
do
./xmain_flops_lapack_orgqr.exe -m $m -n $n -k $k -nb $nb
done
done
done
done
