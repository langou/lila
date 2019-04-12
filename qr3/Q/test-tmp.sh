#!/bin/bash
for i in 32 64 128 256 512 1024 2048 4092
#for i in {1..300}
do
./xmain_flops_mLUinA.exe -n $i
done
