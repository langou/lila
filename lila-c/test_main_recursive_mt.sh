printf "\n\n This is testing mt from 1:10 with ugly everything \n\n"

printf "lila_dgeqrf_w03_mt_l \n"
for i in {1..10}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 0 -panel 0
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {1..10}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 2 -panel 0
done

printf "\n lila_dgeqr2_w03_l \n"
for i in {1..10}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 0
done

printf "\n lila_dgeqr2_w03_3 \n"
for i in {1..10}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 1
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {1..10}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 2
done

#################################################################################################
printf "\n\n This is testing mt from 770:780 - note: mt is larger than $n$ \n\n"

printf "lila_dgeqrf_w03_mt_l \n"
for i in {770..780}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 0 -panel 0
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {770..780}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 2 -panel 0
done

printf "\n lila_dgeqr2_w03_l \n"
for i in {770..780}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 0
done

printf "\n lila_dgeqr2_w03_3 \n"
for i in {770..780}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 1
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {770..780}
do
./xmain.exe -m 1413 -n 692 -mt $i -ii 97  -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 2
done


#################################################################################################
printf "\n\n This is testing (n == m) mt from 1:10 \n\n"

printf "lila_dgeqrf_w03_mt_l \n"
for i in {1..10}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 0 -panel 0
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {1..10}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 2 -panel 0
done

printf "\n lila_dgeqr2_w03_l \n"
for i in {1..10}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 0
done

printf "\n lila_dgeqr2_w03_3 \n"
for i in {1..10}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 1
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {1..10}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 2
done


#################################################################################################
printf "\n\n This is testing (n == m) mt from 770:780 - note: mt is larger than $n$ \n\n"

printf "lila_dgeqrf_w03_mt_l \n"
for i in {770..780}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 0 -panel 0
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {770..780}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 2 -panel 0
done

printf "\n lila_dgeqr2_w03_l \n"
for i in {770..780}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 0
done

printf "\n lila_dgeqr2_w03_3 \n"
for i in {770..780}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 1
done

printf "\n lila_dgeqrf_w03_mt_hr \n"
for i in {770..780}
do
./xmain.exe -m 1213 -n 1200 -ii 13 -mt $i -lda 1902 -ldq 1809 -nx 9 -mode r -leaf 1 -panel 2
done


