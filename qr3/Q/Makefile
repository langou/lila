
include ../make.inc

all: xmain_dorgqr.exe xmain_flops_N2T.exe xmain_flops_V2N.exe xmain_flops_mLUinA.exe xmain_flops_org2r.exe xmain_flops_ApUBTinA.exe xmain_flops_ULTinU.exe xmain_flops_VT2Q.exe xmain_flops_lapack_orgqr.exe

xmain_dorgqr.exe: main_dorgqr.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_larfb.o flops_N2T.o flops_VT2Q.o 
	$(LOADER) -o $@ main_dorgqr.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o our_dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_larfb.o flops_N2T.o flops_VT2Q.o  $(LIB)

xmain_flops_N2T.exe: main_flops_N2T.o flops_trmm.o flops_N2T.o flops_N2T_check.o
	$(LOADER) -o $@ main_flops_N2T.o flops_trmm.o flops_N2T.o flops_N2T_check.o

xmain_flops_V2N.exe: main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o flops_V2N_check.o
	$(LOADER) -o $@ main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o flops_V2N_check.o

xmain_flops_VT2Q.exe: main_flops_VT2Q.o flops_trmm.o flops_ULTinU.o flops_mLUinA.o flops_VT2Q.o flops_VT2Q_check.o
	$(LOADER) -o $@ main_flops_VT2Q.o flops_trmm.o flops_ULTinU.o flops_mLUinA.o flops_VT2Q.o flops_VT2Q_check.o

xmain_flops_mLUinA.exe: main_flops_mLUinA.o flops_trmm.o flops_gemm.o flops_mLUinA.o flops_mLUinA_check.o
	$(LOADER) -o $@ main_flops_mLUinA.o flops_trmm.o flops_gemm.o flops_mLUinA.o flops_mLUinA_check.o

xmain_flops_org2r.exe: main_flops_org2r.o flops_org2r.o
	$(LOADER) -o $@ main_flops_org2r.o flops_org2r.o

xmain_flops_ApUBTinA.exe: main_flops_ApUBTinA.o flops_gemm.o flops_ApUBTinA.o flops_ApUBTinA_check.o
	$(LOADER) -o $@  main_flops_ApUBTinA.o flops_gemm.o flops_ApUBTinA.o flops_ApUBTinA_check.o

xmain_flops_ULTinU.exe: main_flops_ULTinU.o flops_ULTinU.o flops_ULTinU_check.o flops_ApUBTinA.o flops_trmm.o
	$(LOADER) -o $@  main_flops_ULTinU.o flops_ULTinU.o flops_ULTinU_check.o flops_ApUBTinA.o flops_trmm.o

xmain_flops_lapack_orgqr.exe: main_flops_lapack_orgqr.o flops_larft.o flops_larfb.o flops_org2r.o flops_lapack_orgqr.o flops_lapack_orgqr_check.o
	$(LOADER) -o $@  main_flops_lapack_orgqr.o flops_larft.o flops_larfb.o flops_org2r.o flops_lapack_orgqr.o flops_lapack_orgqr_check.o

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe
