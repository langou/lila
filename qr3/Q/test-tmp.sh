#!/bin/bash
for i in {11..20}
do
for j in {1..10}
do
for k in {1..10}
do
./xmain_flops_lapack_orgqr.exe -m $i -n $j -k $k -nb 1
done
done
done
