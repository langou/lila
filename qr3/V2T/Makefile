
include ../make.inc

all: xmain_V2T_lapack.exe 
all: xmain_V2T_qr3.exe 
all: xmain_V2T_qr3_ISW.exe 
all: xmain_V2T_qr3_UT.exe 
all: xmain_V2T_qr3_UT_ISW.exe 
all: xmain_V2T_V2N.exe 
all: xmain_V2T_V2N_N2T.exe 
all: xmain_V2T_V2N_TRI.exe 
all: xmain_V2T_qr3_V2T.exe 
all: xmain_V2T_qr3_V2T_UT.exe 

xmain_V2T_lapack.exe: main_V2T_lapack.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o
	$(LOADER) -o $@ main_V2T_lapack.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o $(LIB)

xmain_V2T_qr3.exe: main_V2T_qr3.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft.o
	$(LOADER) -o $@ main_V2T_qr3.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft.o $(LIB)

xmain_V2T_qr3_ISW.exe: main_V2T_qr3_ISW.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft.o qr3_larft_ISW.o
	$(LOADER) -o $@ main_V2T_qr3_ISW.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft.o qr3_larft_ISW.o $(LIB)

xmain_V2T_qr3_UT.exe: main_V2T_qr3_UT.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft_UT.o qr3_dorgqr_UT.o
	$(LOADER) -o $@ main_V2T_qr3_UT.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft_UT.o qr3_dorgqr_UT.o $(LIB)

xmain_V2T_qr3_UT_ISW.exe: main_V2T_qr3_UT_ISW.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft_UT.o qr3_larft_UT_ISW.o qr3_dorgqr_UT.o
	$(LOADER) -o $@ main_V2T_qr3_UT_ISW.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_larft_UT.o qr3_larft_UT_ISW.o qr3_dorgqr_UT.o $(LIB)

xmain_V2T_V2N.exe: main_V2T_V2N.o qr3_dorgqr_UT.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o dV2N.o dVS2Q.o UinvLTinU.o mLUinA.o
	$(LOADER) -o $@ main_V2T_V2N.o qr3_dorgqr_UT.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o dV2N.o dVS2Q.o UinvLTinU.o mLUinA.o $(LIB)

xmain_V2T_V2N_N2T.exe: main_V2T_V2N_N2T.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o dN2T.o dV2N.o
	$(LOADER) -o $@ main_V2T_V2N_N2T.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o dN2T.o dV2N.o $(LIB)

xmain_V2T_V2N_TRI.exe: main_V2T_V2N_TRI.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o dV2N.o
	$(LOADER) -o $@ main_V2T_V2N_TRI.o qr3_dorgqr.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o dV2N.o $(LIB)

xmain_V2T_qr3_V2T.exe: main_V2T_qr3_V2T.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_dorgqr.o qr3_larft_ISW_V2T.o dN2T.o dV2N.o 
	$(LOADER) -o $@ main_V2T_qr3_V2T.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_dorgqr.o qr3_larft_ISW_V2T.o dN2T.o dV2N.o $(LIB)

xmain_V2T_qr3_V2T_UT.exe: main_V2T_qr3_V2T_UT.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_dorgqr_UT.o qr3_larft_ISW_V2T_UT.o dN2T.o dV2N.o 
	$(LOADER) -o $@ main_V2T_qr3_V2T_UT.o V2T_test_qq_orth_1.o V2T_test_qr_repres_1.o qr3_dorgqr_UT.o qr3_larft_ISW_V2T_UT.o dN2T.o dV2N.o $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe
