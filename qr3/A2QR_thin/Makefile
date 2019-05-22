
include ../make.inc

#all: xmain_A2QR_lapack.exe 
#all: xmain_A2QR_qr3.exe  
#all: xmain_A2QR_qr3_UT.exe  
#all: xmain_A2QR_qr3_UT_ISW.exe  
#all: xmain_A2QR_qr3_ISW.exe  
#all: xmain_A2QR_qr3_VT2Q.exe  
#all: xmain_A2QR_qr3_VT2Q_ISW.exe  
#all: xmain_A2QR_qr3_ISW_VT2Q_ISW.exe  
#all: xmain_A2QR_qr3_UT_VS2Q.exe  
#all: xmain_A2QR_qr3_UT_VS2Q_ISW.exe  

all: xmain_A2QR_lapack.exe xmain_A2QR_qr3.exe xmain_A2QR_qr3_ISW.exe xmain_A2QR_qr3_VT2Q.exe xmain_A2QR_qr3_VT2Q_ISW.exe xmain_A2QR_qr3_UT.exe xmain_A2QR_qr3_ISW_VT2Q_ISW.exe xmain_A2QR_qr3_UT_ISW.exe xmain_A2QR_qr3_UT_VS2Q.exe xmain_A2QR_qr3_UT_VS2Q_ISW.exe  

xmain_A2QR_lapack.exe: main_A2QR_lapack.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_lapack.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3.exe: main_A2QR_qr3.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_ISW.exe: main_A2QR_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_UT.exe: main_A2QR_qr3_UT.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_UT.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_UT_ISW.exe: main_A2QR_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_UT_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_VT2Q.exe: main_A2QR_qr3_VT2Q.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_VT2Q.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_ISW_VT2Q_ISW.exe: main_A2QR_qr3_ISW_VT2Q_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_ISW_VT2Q_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_VT2Q_ISW.exe: main_A2QR_qr3_VT2Q_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_VT2Q_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_UT_VS2Q.exe: main_A2QR_qr3_UT_VS2Q.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_UT_VS2Q.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_UT_VS2Q_ISW.exe: main_A2QR_qr3_UT_VS2Q_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_UT_VS2Q_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

../src/libqr2.a: ../src/qr2_dorgqr_VT2Q.o    \
	../src/qr2_dorgqr_VS2Q.o             \
	../src/qr2_dgeqr3R.o                 \
	../src/qr2_dgeqr3R_UT.o              \
	../src/qr2_dgeqr3R_ISW.o             \
	../src/qr2_dgeqr3R_UT_ISW.o          \
	../src/qr2_dorgqr.o                  \
	../src/qr2_dorgqr_UT.o               \
	../src/qr2_dVT2Q.o                   \
	../src/qr2_dVS2Q.o                   \
	../src/qr2_ULTinU.o                  \
	../src/qr2_mLUinA.o              
	make -C ../src/

../check/libcheck.a: ../check/check_qq_orth.o \
	../check/check_qr_repres.o
	make -C ../check/

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o
