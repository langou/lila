#!/bin/bash
for i in {1..1000}
do
./xmain.exe -mode recursive -w 3 -m 1000 -n 1000 -mt $i -lda 1000
done
