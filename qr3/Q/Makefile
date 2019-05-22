
include ../make.inc

#all: xmain_A2QR_ourlapack.exe xmain_A2QR_lapack.exe xmain_A2QR_NminusKblock.exe xmain_A2QRTV.exe xmain_study_of_larft.exe xmain_A2QR_ISW.exe xmain_A2QR_UT.exe xmain_A2QR_ISW_UT.exe xmain_A2Q1Q2R.exe

#all: xmain_A2QR_lapack.exe
#all: xmain_A2QR_ourlapack.exe
#all: xmain_A2QR_NminusKblock.exe
#all: xmain_A2QRTV.exe
#all: xmain_A2QR_ISW.exe 
#all: xmain_A2QR_UT.exe
#all: xmain_A2QR_ISW_UT.exe 
all: xmain_A2Q1Q2R.exe 
#all: xmain_study_of_larft.exe

xmain_A2QR_lapack.exe: main_A2QR_lapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o 
	$(LOADER) -o $@ main_A2QR_lapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o $(LIB)

xmain_A2QR_NminusKblock.exe: main_A2QR_NminusKblock.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dgeqrf.o dorgqr_after.o dgeqr3.o  dgeqr3.o
	$(LOADER) -o $@ main_A2QR_NminusKblock.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dgeqrf.o dorgqr_after.o dgeqr3.o $(LIB)

xmain_A2QR_ourlapack.exe: main_A2QR_ourlapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dgeqrf.o dorgqr_after.o qr3_dorgqr_level1.o dgeqr3R.o qr3_dorgqr_level1_UT.o dVS2Q.o UinvLTinU.o UinvLTinU_cheat.o qr3_dorgqr.o qr3_larft.o
	$(LOADER) -o $@ main_A2QR_ourlapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dgeqrf.o dorgqr_after.o qr3_dorgqr_level1.o dgeqr3R.o qr3_dorgqr_level1_UT.o dVS2Q.o UinvLTinU.o UinvLTinU_cheat.o qr3_dorgqr.o qr3_larft.o $(LIB)

xmain_study_of_larft.exe: main_study_of_larft.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dgeqrf.o dorgqr_after.o qr3_dorgqr_level1.o dgeqr3.o qr3_dorgqr_level1_UT.o dVS2Q.o UinvLTinU.o UinvLTinU_cheat.o qr3_dorgqr.o qr3_larft.o
	$(LOADER) -o $@ main_study_of_larft.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dgeqrf.o dorgqr_after.o qr3_dorgqr_level1.o dgeqr3.o qr3_dorgqr_level1_UT.o dVS2Q.o UinvLTinU.o UinvLTinU_cheat.o qr3_dorgqr.o qr3_larft.o $(LIB)


xmain_A2QR_UT.exe: main_A2QR_UT.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o qr3_dorgqr_level1_UT.o qr3_dorgqr_level1.o dN2T.o dV2N.o dVS2Q.o UinvLTinU.o mLUinA.o 
	$(LOADER) -o $@ main_A2QR_UT.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o qr3_dorgqr_level1_UT.o qr3_dorgqr_level1.o dN2T.o dV2N.o dVS2Q.o UinvLTinU.o mLUinA.o $(LIB)

xmain_A2QR_ISW.exe: main_A2QR_ISW.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dgeqr3.o dgeqr3_ISW.o dVT2Q.o qr3_dorgqr.o ULTinU.o mLUinA.o ApUBTinA.o
	$(LOADER) -o $@ main_A2QR_ISW.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dgeqr3.o dgeqr3_ISW.o dVT2Q.o qr3_dorgqr.o ULTinU.o mLUinA.o ApUBTinA.o $(LIB)

xmain_A2QR_ISW_UT.exe: main_A2QR_ISW_UT.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dgeqr3.o dgeqr3_ISW.o dVT2Q.o qr3_dorgqr.o ULTinU.o mLUinA.o ApUBTinA.o qr3_dorgqr_level1.o dN2T.o dV2N.o qr3_dorgqr_level1_UT.o dVS2Q.o UinvLTinU.o dgeqr3R.o
	$(LOADER) -o $@ main_A2QR_ISW_UT.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dgeqr3.o dgeqr3_ISW.o dVT2Q.o qr3_dorgqr.o ULTinU.o mLUinA.o ApUBTinA.o qr3_dorgqr_level1.o dN2T.o dV2N.o qr3_dorgqr_level1_UT.o dVS2Q.o UinvLTinU.o dgeqr3R.o $(LIB)

xmain_A2QRTV.exe: main_A2QRTV.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o qr3_dA2QRTV_fake.o dgeqr3.o dgeqr3_ISW.o qr3_dorgqr.o qr3_dA2QRTV.o our_dlarfb_lnfc.o qr3_larft.o
	$(LOADER) -o $@ main_A2QRTV.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o qr3_dA2QRTV_fake.o dgeqr3.o dgeqr3_ISW.o qr3_dorgqr.o qr3_dA2QRTV.o our_dlarfb_lnfc.o qr3_larft.o $(LIB)

xmain_A2Q1Q2R.exe: main_A2Q1Q2R.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dgeqr3.o dgeqr3R.o dgeqr3_ISW.o qr3_dorgqr.o qr3_null_dorgqr.o dN2T.o dV2N.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o dorgqr_after.o our_dorg2r_Q2.o our_dorgqr_Q2.o wrapper_dlarf.o our_dorg2r.o
	$(LOADER) -o $@ main_A2Q1Q2R.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dgeqr3.o dgeqr3R.o dgeqr3_ISW.o qr3_dorgqr.o qr3_null_dorgqr.o dN2T.o dV2N.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o dorgqr_after.o our_dorg2r_Q2.o our_dorgqr_Q2.o wrapper_dlarf.o our_dorg2r.o $(LIB)

our_dorgqr_fortran.o: our_dorgqr_fortran.f
	gfortran -c our_dorgqr_fortran.f 

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe
