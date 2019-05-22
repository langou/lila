
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

xmain_V2T_lapack.exe: main_V2T_lapack.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_lapack.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3.exe: main_V2T_qr3.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_qr3.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_ISW.exe: main_V2T_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_UT.exe: main_V2T_qr3_UT.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_qr3_UT.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_UT_ISW.exe: main_V2T_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_V2N.exe: main_V2T_V2N.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_V2N.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_V2N_N2T.exe: main_V2T_V2N_N2T.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_V2N_N2T.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_V2N_TRI.exe: main_V2T_V2N_TRI.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_V2N_TRI.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_V2T.exe: main_V2T_qr3_V2T.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_qr3_V2T.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_V2T_qr3_V2T_UT.exe: main_V2T_qr3_V2T_UT.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a
	$(LOADER) -o $@ main_V2T_qr3_V2T_UT.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

../src/libqr2.a: ../src/qr2_dorgqr.o  \
	../src/qr2_dorgqr_UT.o        \
	../src/qr2_UinvLTinU.o        \
	../src/qr2_dVS2Q.o            \
	../src/qr2_dV2N.o             \
	../src/qr2_dN2T.o             \
	../src/qr2_larft.o            \
	../src/qr2_larft_UT.o         \
	../src/qr2_larft_ISW.o        \
	../src/qr2_larft_ISW_V2T.o    \
	../src/qr2_larft_ISW_V2T_UT.o \
	../src/qr2_larft_UT_ISW.o   
	make -C ../src/

../check/libcheck.a: ../check/check_qr_repres.o \
	../check/check_qq_orth.o
	make -C ../check/

../flops/libflops.a: ../flops/flops_larft.o      
	make -C ../flops/

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o
