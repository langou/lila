
include ../make.inc

all: xmain_A2QR_ourlapack.exe xmain_A2QR_lapack.exe xmain_flops_N2T.exe xmain_flops_V2N.exe xmain_flops_mLUinA.exe xmain_flops_org2r.exe xmain_flops_ApUBTinA.exe xmain_flops_ULTinU.exe xmain_flops_VT2Q.exe xmain_flops_lapack_orgqr.exe xmain_flops_lapack_larfb.exe xmain_flops_lapack_geqrf.exe

xmain_A2QR_lapack.exe: main_A2QR_lapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o
	$(LOADER) -o $@ main_A2QR_lapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o  $(LIB)

xmain_A2QR_ourlapack.exe: main_A2QR_ourlapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o
	$(LOADER) -o $@ main_A2QR_ourlapack.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_N2T.o flops_VT2Q.o our_dlarfb_lnfc.o our_dorgqr_fortran.o our_dgeqrf.o $(LIB)

xmain_flops_N2T.exe: main_flops_N2T.o flops_trmm.o flops_N2T.o flops_N2T_check.o
	$(LOADER) -o $@ main_flops_N2T.o flops_trmm.o flops_N2T.o flops_N2T_check.o

xmain_flops_V2N.exe: main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o flops_V2N_check.o
	$(LOADER) -o $@ main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o flops_V2N_check.o

xmain_flops_VT2Q.exe: main_flops_VT2Q.o flops_trmm.o flops_ULTinU.o flops_mLUinA.o flops_VT2Q.o flops_VT2Q_check.o
	$(LOADER) -o $@ main_flops_VT2Q.o flops_trmm.o flops_ULTinU.o flops_mLUinA.o flops_VT2Q.o flops_VT2Q_check.o

xmain_flops_mLUinA.exe: main_flops_mLUinA.o flops_trmm.o flops_gemm.o flops_mLUinA.o flops_mLUinA_check.o
	$(LOADER) -o $@ main_flops_mLUinA.o flops_trmm.o flops_gemm.o flops_mLUinA.o flops_mLUinA_check.o

xmain_flops_org2r.exe: main_flops_org2r.o flops_org2r.o flops_org2r_check.o flops_lapack_larf.o flops_lapack_larfb.o
	$(LOADER) -o $@ main_flops_org2r.o flops_org2r.o flops_org2r_check.o flops_lapack_larf.o flops_lapack_larfb.o

xmain_flops_ApUBTinA.exe: main_flops_ApUBTinA.o flops_gemm.o flops_ApUBTinA.o flops_ApUBTinA_check.o
	$(LOADER) -o $@  main_flops_ApUBTinA.o flops_gemm.o flops_ApUBTinA.o flops_ApUBTinA_check.o

xmain_flops_ULTinU.exe: main_flops_ULTinU.o flops_ULTinU.o flops_ULTinU_check.o flops_ApUBTinA.o flops_trmm.o
	$(LOADER) -o $@  main_flops_ULTinU.o flops_ULTinU.o flops_ULTinU_check.o flops_ApUBTinA.o flops_trmm.o

xmain_flops_lapack_orgqr.exe: main_flops_lapack_orgqr.o flops_larft.o flops_lapack_larfb.o flops_org2r.o flops_lapack_orgqr.o flops_lapack_orgqr_check.o flops_lapack_orgqr_from_last_org2r.o flops_lapack_orgqr_from_first_org2r.o flops_lapack_orgqr_from_larft.o  flops_org2r_n1.o
	$(LOADER) -o $@  main_flops_lapack_orgqr.o flops_larft.o flops_lapack_larfb.o flops_org2r.o flops_lapack_orgqr.o flops_lapack_orgqr_check.o flops_lapack_orgqr_from_last_org2r.o flops_lapack_orgqr_from_first_org2r.o  flops_lapack_orgqr_from_larft.o flops_org2r_n1.o flops_lapack_larf.o

xmain_flops_lapack_larfb.exe: main_flops_lapack_larfb.o flops_gemm.o flops_trmm.o flops_lapack_larfb.o flops_lapack_larfb_check.o
	$(LOADER) -o $@  main_flops_lapack_larfb.o flops_gemm.o flops_trmm.o flops_lapack_larfb.o flops_lapack_larfb_check.o

xmain_flops_lapack_geqrf.exe: main_flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_geqrf.o flops_lapack_geqrf_check.o flops_lapack_larf.o flops_lapack_geqr2.o flops_larft.o flops_lapack_larfg.o
	$(LOADER) -o $@  main_flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_geqrf.o flops_lapack_geqrf_check.o flops_lapack_larf.o flops_lapack_geqr2.o flops_larft.o flops_lapack_larfg.o

our_dorgqr_fortran.o: our_dorgqr_fortran.f
	gfortran -c our_dorgqr_fortran.f 

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe
