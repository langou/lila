
include ../make.inc

#all: xmain_A2VR_lapack.exe 
#all: xmain_A2VR_qr3.exe 
#all: xmain_A2VR_qr3_ISW.exe 
#all: xmain_A2VR_qr3_UT.exe 
#all: xmain_A2VR_qr3_UT_ISW.exe 

all: xmain_A2VR_lapack.exe xmain_A2VR_qr3.exe xmain_A2VR_qr3_ISW.exe xmain_A2VR_qr3_UT.exe xmain_A2VR_qr3_UT_ISW.exe 

xmain_A2VR_lapack.exe: main_A2VR_lapack.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2VR_lapack.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2VR_qr3.exe: main_A2VR_qr3.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2VR_qr3.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2VR_qr3_ISW.exe: main_A2VR_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2VR_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2VR_qr3_UT.exe: main_A2VR_qr3_UT.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2VR_qr3_UT.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2VR_qr3_UT_ISW.exe: main_A2VR_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2VR_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

../src/libqr2.a: ../src/qr2_dgeqr3R.o \
	../src/qr2_dgeqr3R_ISW.o      \
	../src/qr2_dgeqr3R_UT.o       \
	../src/qr2_dgeqr3R_UT_ISW.o   \
	../src/qr2_dorgqr.o           \
	../src/qr2_dV2tau.o   
	make -C ../src/

../check/libcheck.a: ../check/check_qr_repres.o \
	../check/check_qq_orth.o
	make -C ../check/

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o
