#!/bin/sh
for i in {1..402}
do
./xmain_unittest__lila_dgeqrf_w03_l -m 701 -n 402 -mt $i -nb 77 -ii 10 -lda 902 -ldq 809
done
