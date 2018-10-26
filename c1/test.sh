#!/bin/sh
./xmain.exe -mode recursive -w 2 -m 100 -n 49 -mt 7
./xmain.exe -mode recursive -w 3 -m 100 -n 49 -mt 7
./xmain.exe -mode levelx    -w 2 -m 100 -n 49 -mt 7 -n_lvl 1 3
./xmain.exe -mode levelx    -w 3 -m 100 -n 49 -mt 7 -n_lvl 1 3

./xmain.exe -mode recursive -w 2 -m 400 -n 400 -mt 15
./xmain.exe -mode recursive -w 3 -m 400 -n 400 -mt 15
./xmain.exe -mode levelx    -w 2 -m 400 -n 400 -mt 15 -n_lvl 2 116 33
./xmain.exe -mode levelx    -w 3 -m 400 -n 400 -mt 15 -n_lvl 2 116 33
./xmain.exe -mode levelx    -w 2 -m 400 -n 400 -mt 15 -n_lvl 3 73 11 2
./xmain.exe -mode levelx    -w 3 -m 400 -n 400 -mt 15 -n_lvl 3 73 11 2
