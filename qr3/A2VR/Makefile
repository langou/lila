
include ../make.inc

all: xmain_A2VR_lapack.exe 
all: xmain_A2VR_qr3.exe 
all: xmain_A2VR_qr3_ISW.exe 
all: xmain_A2VR_qr3_UT.exe 
all: xmain_A2VR_qr3_UT_ISW.exe 

xmain_A2VR_lapack.exe: main_A2VR_lapack.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o
	$(LOADER) -o $@ main_A2VR_lapack.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o $(LIB)

xmain_A2VR_qr3.exe: main_A2VR_qr3.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R.o dV2tau.o
	$(LOADER) -o $@ main_A2VR_qr3.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R.o dV2tau.o $(LIB)

xmain_A2VR_qr3_ISW.exe: main_A2VR_qr3_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R.o dgeqr3R_ISW.o dV2tau.o
	$(LOADER) -o $@ main_A2VR_qr3_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R.o dgeqr3R_ISW.o dV2tau.o $(LIB)

xmain_A2VR_qr3_UT.exe: main_A2VR_qr3_UT.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dV2tau.o
	$(LOADER) -o $@ main_A2VR_qr3_UT.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dV2tau.o $(LIB)

xmain_A2VR_qr3_UT_ISW.exe: main_A2VR_qr3_UT_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dgeqr3R_UT_ISW.o dV2tau.o 
	$(LOADER) -o $@ main_A2VR_qr3_UT_ISW.o qr3_dorgqr.o A2VR_test_qq_orth_1.o A2VR_test_qr_repres_1.o dgeqr3R_UT.o dgeqr3R_UT_ISW.o dV2tau.o $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe
