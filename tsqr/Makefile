
all: xmain_test_tsqr.exe 

xmain_test_tsqr.exe: main_test_tsqr.o dlamtsqr.o dorhrhack.o dlatsqr.o
	gfortran -o $@ main_test_tsqr.o  dlamtsqr.o dorhrhack.o dlatsqr.o  /Users/langou/Desktop/repositories/lapack.git/liblapacke.a /Users/langou/Desktop/repositories/lapack.git/liblapack.a /Users/langou/Desktop/repositories/lapack.git/libcblas.a /Users/langou/Desktop/repositories/lapack.git/librefblas.a

dlamtsqr.o: dlamtsqr.f
	gfortran -c dlamtsqr.f

dlatsqr.o: dlatsqr.f
	gfortran -c dlatsqr.f

dorhrhack.o: dorhrhack.f
	gfortran -c dorhrhack.f

.c.o:
	gcc -I/Users/langou/Desktop/repositories/lapack.git/LAPACKE/include -I/Users/langou/Desktop/repositories/lapack.git/CBLAS/include -c -o $@ $<

clean:
	rm -f *o *exe
