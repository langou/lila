#!/bin/bash
for i in {1..10}
do
./xmain_flops_ULTinU.exe -n $i
done

for i in {1..10}
do
for j in {1..10}
do
./xmain_flops_ApUBTinA.exe -m $i -n $j
done
done

for i in {1..10}
do
./xmain_flops_V2N.exe -n $i
done

for i in {1..10}
do
./xmain_flops_mLUinA.exe -n $i
done

for i in {1..10}
do
./xmain_flops_N2T.exe -n $i
done





#	dN2T( ib, tauk, Akk, lda );
#	dVT2Q( ml, ib, Akk, lda );

