
include ../make.inc

all: xmain_test__lapack_Q2.exe xmain_test__lapack_mod_Q2.exe

xmain_test__lapack_Q2.exe: main_test__lapack_Q2.o 
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_test__lapack_Q2.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_test__lapack_mod_Q2.exe: main_test__lapack_mod_Q2.o 
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_test__lapack_mod_Q2.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o ../flops/*a ../flops/*o
