#!/bin/bash
#for i in 32 64 128 256 512 1024 2048 4092
for i in {1..300}
do
./a.out -n $i
done
