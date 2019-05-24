
include ../make.inc

#all: xmain_V2T_lapack.exe 
#all: xmain_V2T_qr3.exe 
#all: xmain_V2T_qr3_ISW.exe 
#all: xmain_V2T_qr3_UT.exe 
#all: xmain_V2T_qr3_UT_ISW.exe 
#all: xmain_V2T_V2N.exe 
#all: xmain_V2T_V2N_N2T.exe 
#all: xmain_V2T_V2N_TRI.exe 
#all: xmain_V2T_qr3_V2T.exe 
#all: xmain_V2T_qr3_V2T_UT.exe 

all: xmain_V2T_lapack.exe xmain_V2T_qr3.exe xmain_V2T_qr3_ISW.exe xmain_V2T_qr3_UT.exe xmain_V2T_qr3_UT_ISW.exe xmain_V2T_V2N.exe xmain_V2T_V2N_N2T.exe xmain_V2T_V2N_TRI.exe xmain_V2T_qr3_V2T.exe xmain_V2T_qr3_V2T_UT.exe 

xmain_V2T_lapack.exe: main_V2T_lapack.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_lapack.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3.exe: main_V2T_qr3.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_qr3.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_ISW.exe: main_V2T_qr3_ISW.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_UT.exe: main_V2T_qr3_UT.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_qr3_UT.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_UT_ISW.exe: main_V2T_qr3_UT_ISW.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_V2N.exe: main_V2T_V2N.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_V2N.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_V2N_N2T.exe: main_V2T_V2N_N2T.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_V2N_N2T.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_V2N_TRI.exe: main_V2T_V2N_TRI.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_V2N_TRI.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_V2T.exe: main_V2T_qr3_V2T.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_qr3_V2T.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_V2T_UT.exe: main_V2T_qr3_V2T_UT.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_V2T_qr3_V2T_UT.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o
