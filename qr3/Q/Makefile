
include ../make.inc

#all: xmain_A2QR_ourlapack.exe xmain_A2QR_lapack.exe xmain_flops_N2T.exe xmain_flops_V2N.exe xmain_flops_mLUinA.exe xmain_flops_ApUBTinA.exe xmain_flops_ULTinU.exe xmain_flops_VT2Q.exe xmain_flops_lapack_orgqr.exe xmain_flops_lapack_larfb.exe xmain_flops_lapack_geqrf.exe xmain_A2QR_NminusKblock.exe xmain_flops_dorgqr_after.exe xmain_A2QR_ISW.exe xmain_A2QRTV.exe xmain_flops_lapack_geqr2.exe

all: xmain_flops_lapack_orgqr.exe xmain_flops_lapack_org2r.exe xmain_flops_legacy_lapack_org2r.exe

xmain_A2QR_lapack.exe: main_A2QR_lapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o
	$(LOADER) -o $@ main_A2QR_lapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o  $(LIB)

xmain_A2QR_ourlapack.exe: main_A2QR_ourlapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o dorgqr_after.o
	$(LOADER) -o $@ main_A2QR_ourlapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o dorgqr_after.o $(LIB)

xmain_A2QR_NminusKblock.exe: main_A2QR_NminusKblock.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o dorgqr_after.o
	$(LOADER) -o $@ main_A2QR_NminusKblock.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o dorgqr_after.o $(LIB)

xmain_A2QR_ISW.exe: main_A2QR_ISW.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o dorgqr_after.o dgeqr3_right.o dgeqr3_left.o qr3_larft.o qr3_dorgqr.o
	$(LOADER) -o $@ main_A2QR_ISW.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o dorgqr_after.o dgeqr3_right.o dgeqr3_left.o qr3_larft.o qr3_dorgqr.o $(LIB)

xmain_A2QRTV.exe: main_A2QRTV.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o qr3_dA2QRTV_fake.o dgeqr3_left.o dgeqr3_right.o qr3_dorgqr.o qr3_dA2QRTV.o our_dlarfb_lnfc.o qr3_larft.o
	$(LOADER) -o $@ main_A2QRTV.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o qr3_dA2QRTV_fake.o dgeqr3_left.o dgeqr3_right.o qr3_dorgqr.o qr3_dA2QRTV.o our_dlarfb_lnfc.o qr3_larft.o $(LIB)

xmain_flops_N2T.exe: main_flops_N2T.o flops_trmm.o flops_N2T.o flops_N2T_check.o
	$(LOADER) -o $@ main_flops_N2T.o flops_trmm.o flops_N2T.o flops_N2T_check.o

xmain_flops_V2N.exe: main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o flops_V2N_check.o
	$(LOADER) -o $@ main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o flops_V2N_check.o

xmain_flops_VT2Q.exe: main_flops_VT2Q.o flops_trmm.o flops_ULTinU.o flops_mLUinA.o flops_VT2Q.o flops_VT2Q_check.o
	$(LOADER) -o $@ main_flops_VT2Q.o flops_trmm.o flops_ULTinU.o flops_mLUinA.o flops_VT2Q.o flops_VT2Q_check.o

xmain_flops_mLUinA.exe: main_flops_mLUinA.o flops_trmm.o flops_gemm.o flops_mLUinA.o flops_mLUinA_check.o
	$(LOADER) -o $@ main_flops_mLUinA.o flops_trmm.o flops_gemm.o flops_mLUinA.o flops_mLUinA_check.o

xmain_flops_ApUBTinA.exe: main_flops_ApUBTinA.o flops_gemm.o flops_ApUBTinA.o flops_ApUBTinA_check.o
	$(LOADER) -o $@  main_flops_ApUBTinA.o flops_gemm.o flops_ApUBTinA.o flops_ApUBTinA_check.o

xmain_flops_ULTinU.exe: main_flops_ULTinU.o flops_ULTinU.o flops_ULTinU_check.o flops_ApUBTinA.o flops_trmm.o
	$(LOADER) -o $@  main_flops_ULTinU.o flops_ULTinU.o flops_ULTinU_check.o flops_ApUBTinA.o flops_trmm.o

xmain_flops_lapack_larfb.exe: main_flops_lapack_larfb.o flops_gemm.o flops_trmm.o flops_lapack_larfb.o flops_lapack_larfb_check.o
	$(LOADER) -o $@  main_flops_lapack_larfb.o flops_gemm.o flops_trmm.o flops_lapack_larfb.o flops_lapack_larfb_check.o

xmain_flops_lapack_geqrf.exe: main_flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_geqrf.o flops_lapack_geqrf_check.o flops_lapack_larf.o flops_lapack_geqr2.o flops_larft.o flops_lapack_larfg.o flops_lapack_geqr2_check.o
	$(LOADER) -o $@  main_flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_geqrf.o flops_lapack_geqrf_check.o flops_lapack_larf.o flops_lapack_geqr2.o flops_larft.o flops_lapack_larfg.o flops_lapack_geqr2_check.o

xmain_flops_lapack_geqr2.exe: main_flops_lapack_geqr2.o flops_lapack_geqr2.o flops_lapack_geqr2_check.o flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_larfg.o flops_larft.o flops_lapack_larf.o flops_lapack_geqrf_check.o
	$(LOADER) -o $@  main_flops_lapack_geqr2.o flops_lapack_geqr2.o flops_lapack_geqr2_check.o flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_larfg.o flops_larft.o flops_lapack_larf.o flops_lapack_geqrf_check.o

xmain_flops_lapack_orgqr.exe: main_flops_lapack_orgqr.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o flops_lapack_larfb.o flops_lapack_org2r.o
	$(LOADER) -o $@  main_flops_lapack_orgqr.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o  flops_lapack_larfb.o flops_lapack_org2r.o


xmain_flops_lapack_org2r.exe: main_flops_lapack_org2r.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o flops_lapack_larfb.o flops_lapack_org2r.o
	$(LOADER) -o $@  main_flops_lapack_org2r.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o  flops_lapack_larfb.o flops_lapack_org2r.o

xmain_flops_legacy_lapack_org2r.exe: main_flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r_check.o flops_legacy_lapack_larf.o
	$(LOADER) -o $@  main_flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r_check.o flops_legacy_lapack_larf.o

our_dorgqr_fortran.o: our_dorgqr_fortran.f
	gfortran -c our_dorgqr_fortran.f 

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe
