#!/bin/sh
for i in {1..402}
do
./xmain_unittest__lila_dgeqrf_w03_l -m 1213 -n 692 -mt $i -ii 101  -lda 1902 -ldq 1809
done
