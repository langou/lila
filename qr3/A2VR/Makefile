
include ../make.inc

#all: xmain_A2VR_lapack.exe 
#all: xmain_A2VR_qr3.exe 
#all: xmain_A2VR_qr3_ISW.exe 
#all: xmain_A2VR_qr3_UT.exe 
#all: xmain_A2VR_qr3_UT_ISW.exe 

#all: xmain_A2VR_lapack.exe xmain_A2VR_qr3.exe xmain_A2VR_qr3_ISW.exe xmain_A2VR_qr3_UT.exe xmain_A2VR_qr3_UT_ISW.exe 

xmain_A2VR_lapack.exe: main_A2VR_lapack.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o ../src/libqr2.a
	$(LOADER) -o $@ main_A2VR_lapack.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o ../src/libqr2.a $(LIB)

xmain_A2VR_qr3.exe: main_A2VR_qr3.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dV2tau.o ../src/libqr2.a
	$(LOADER) -o $@ main_A2VR_qr3.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dV2tau.o ../src/libqr2.a $(LIB)

xmain_A2VR_qr3_ISW.exe: main_A2VR_qr3_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dV2tau.o ../src/libqr2.a
	$(LOADER) -o $@ main_A2VR_qr3_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dV2tau.o ../src/libqr2.a $(LIB)

xmain_A2VR_qr3_UT.exe: main_A2VR_qr3_UT.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dV2tau.o ../src/libqr2.a
	$(LOADER) -o $@ main_A2VR_qr3_UT.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dV2tau.o ../src/libqr2.a $(LIB)

xmain_A2VR_qr3_UT_ISW.exe: main_A2VR_qr3_UT_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dgeqr3R_UT_ISW.o dV2tau.o ../src/libqr2.a
	$(LOADER) -o $@ main_A2VR_qr3_UT_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dgeqr3R_UT_ISW.o dV2tau.o ../src/libqr2.a $(LIB)

../src/libqr2.a: ../src/qr2_dgeqr3R.o \
	../src/qr2_dgeqr3R_ISW.o
	make -C ../src/

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o
