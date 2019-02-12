#!/bin/sh
#./xmain.exe -mode recursive -w 2 -m 100 -n 49 -mt 7
#./xmain.exe -mode recursive -w 3 -m 100 -n 49 -mt 7
#./xmain.exe -mode levelx    -w 2 -m 100 -n 49 -mt 7 -n_lvl 1 3
#./xmain.exe -mode levelx    -w 3 -m 100 -n 49 -mt 7 -n_lvl 1 3

#./xmain.exe -mode recursive -w 2 -m 400 -n 400 -mt 15
#./xmain.exe -mode recursive -w 3 -m 400 -n 400 -mt 15
#./xmain.exe -mode levelx    -w 2 -m 400 -n 400 -mt 15 -n_lvl 2 116 33
#./xmain.exe -mode levelx    -w 3 -m 400 -n 400 -mt 15 -n_lvl 2 116 33
#./xmain.exe -mode levelx    -w 2 -m 400 -n 400 -mt 15 -n_lvl 3 73 11 2
#./xmain.exe -mode levelx    -w 3 -m 400 -n 400 -mt 15 -n_lvl 3 73 11 2

./xmain.exe -mode recursive -w 3 -m 1213 -n 692 -mt 1 -ii 101  -lda 1902 -ldq 1809
for i in {200..402}
do
./xmain.exe -mode recursive -w 3 -m 1213 -n 692 -mt $i -ii 101  -lda 1902 -ldq 1809
done
./xmain.exe -mode recursive -w 3 -m 1213 -n 692 -mt 1213 -ii 101  -lda 1902 -ldq 1809
for i in {200..303}
do
./xmain.exe -mode recursive -w 3 -m 1213 -n 692 -mt 101 -ii $i  -lda 1902 -ldq 1809
done



./xmain.exe -mode levelx -w 3 -m 1213 -n 692 -mt 101 -ii $i  -lda 1902 -ldq 1809
for i in {200..402}
do
./xmain.exe -mode levelx -w 3 -m 1213 -n 692 -mt $i -ii 101  -lda 1902 -ldq 1809
done
./xmain.exe -mode levelx -w 3 -m 1213 -n 692 -mt $i -ii 101  -lda 1902 -ldq 1809
for i in {200..303}
do
./xmain.exe -mode levelx -w 3 -m 1213 -n 692 -mt 101 -ii $i  -lda 1902 -ldq 1809
done
