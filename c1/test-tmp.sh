#!/bin/bash
#./xmain.exe -mode recursive -w 3 -m 20 -n 19 -mt 3
#./xmain.exe -mode recursive -w 3 -m 100 -n 39 -mt 3
#for i in {1..420}
#do
#./xmain.exe -mode recursive -w 3 -m 1000 -n 400 -mt $i -lda 1201
#./xmain.exe -mode levelx    -w 3 -m 1000 -n 400 -mt $i -lda 1201 -n_lvl 1 10
#./xmain.exe -mode levelx    -w 3 -m 1000 -n 400 -mt $i -lda 1201 -n_lvl 3 103 27 4
#done
for i in {1..300}
do
./xmain.exe -mode recursive -w 3 -m 1000 -n 501 -mt $i
done
