#!/bin/sh
for i in {1..499}
do
./xmain_demo_03b.exe -m 500 -n 499 -mt $i -nb 14
done
