
all: xmain_test_tsqr.exe 

xmain_test_tsqr.exe: main_test_tsqr.o
	gfortran -o $@ main_test_tsqr.o  /Users/langou/Desktop/repositories/lapack.git/liblapacke.a /Users/langou/Desktop/repositories/lapack.git/liblapack.a /Users/langou/Desktop/repositories/lapack.git/libcblas.a /Users/langou/Desktop/repositories/lapack.git/librefblas.a

.c.o:
	gcc -I/Users/langou/Desktop/repositories/lapack.git/LAPACKE/include -I/Users/langou/Desktop/repositories/lapack.git/CBLAS/include -c -o $@ $<

clean:
	rm -f *o *exe
