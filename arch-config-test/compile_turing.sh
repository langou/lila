
echo 'rm -f main_dgemm.o'
rm -f main_dgemm.o

echo 'gcc -c -D_USE_BLIS -I/home/math/langou/opt/blis.git/include/haswell/ main_dgemm.c -o main_dgemm.o'
gcc -c -D_USE_BLIS -I/home/math/langou/opt/blis.git/include/haswell/ main_dgemm.c -o main_dgemm.o

echo 'gcc main_dgemm.o /home/math/langou/opt/blis.git/lib/haswell/libblis.a -lm -lpthread -o xmain_dgemm_blis.exe'
gcc main_dgemm.o /home/math/langou/opt/blis.git/lib/haswell/libblis.a -lm -lpthread -o xmain_dgemm_blis.exe

echo 'rm -f main_dgemm.o'
rm -f main_dgemm.o

echo 'gcc -c -D_USE_MKL -I/home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/include/ main_dgemm.c -o main_dgemm.o'
gcc -c -D_USE_MKL -I/home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/include/ main_dgemm.c -o main_dgemm.o

echo 'gcc main_dgemm.o -Wl,--start-group /home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/libmkl_intel_lp64.a /home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/libmkl_gnu_thread.a /home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -o xmain_dgemm_mkl.exe'
gcc main_dgemm.o -Wl,--start-group /home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/libmkl_intel_lp64.a /home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/libmkl_gnu_thread.a /home/math/langou/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -o xmain_dgemm_mkl.exe

echo 'rm -f main_dgemm.o'
rm -f main_dgemm.o

echo 'gcc -c -I/home/math/langou/opt/openblas.git/ main_dgemm.c -o main_dgemm.o'
gcc -c -I/home/math/langou/opt/openblas.git/ main_dgemm.c -o main_dgemm.o

echo 'gcc main_dgemm.o /home/math/langou/opt/openblas.git/libopenblas.a -lpthread -o xmain_dgemm_openblas.exe'
gcc main_dgemm.o /home/math/langou/opt/openblas.git/libopenblas.a -lpthread -o xmain_dgemm_openblas.exe

echo 'rm -f main_dgemm.o'
rm -f main_dgemm.o

