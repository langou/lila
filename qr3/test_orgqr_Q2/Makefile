
include ../make.inc

all: xmain_test_orgqr_Q2.exe 

xmain_test_orgqr_Q2.exe: main_test_orgqr_Q2.o qr3_aux_dorgq2r.o lapack_our_dorgq2r.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_test_orgqr_Q2.o qr3_aux_dorgq2r.o lapack_our_dorgq2r.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o ../flops/*a ../flops/*o
