#!/bin/bash
for i in {1..50}
do
for j in {1..50}
do
./xmain_flops_ApUBTinA.exe -m $i -n $j
done
done
