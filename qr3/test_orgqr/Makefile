
include ../make.inc

#all: xmain_test_orgqr__lapack.exe 
#all: xmain_test_orgqr__lapack_ref.exe 
all: xmain_test_orgqr__lapack_mod.exe 

xmain_test_orgqr__lapack.exe: main_test_orgqr__lapack.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_test_orgqr__lapack.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_test_orgqr__lapack_ref.exe: main_test_orgqr__lapack_ref.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_test_orgqr__lapack_ref.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_test_orgqr__lapack_mod.exe: main_test_orgqr__lapack_mod.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_test_orgqr__lapack_mod.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o ../flops/*a ../flops/*o
