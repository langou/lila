
include ../make.inc

FLOPS = flops_geqr3.o      \
	flops_geqr3_ISW.o  \

CC = gcc
OPTS = -Wall

all: $(FLOPS)
	ar -cr libflops.a $^
	ranlib libflops.a

#all: xmain_flops_lapack_orgqr.exe xmain_flops_lapack_org2r.exe xmain_flops_lapack_geqr2.exe xmain_flops_lapack_geqrf.exe

#all: xmain_flops_geqr3.exe 

#all: xmain_flops_geqr3_ISW.exe

xmain_flops_geqr3_ISW.exe: main_flops_geqr3_ISW.o flops_geqr3_ISW_check.o flops_gemm.o flops_trmm.o flops_geqr3_check.o flops_lapack_larfg.o flops_geqr3_ISW.o flops_geqr3_ISW_save.o flops_geqr3.o flops_larft.o flops_lapack_geqr2.o flops_geqr3_bef_constructT_check.o flops_geqr3_bef_useT.o flops_geqr3_bef_useT_check.o flops_geqr3_wob_check.o flops_geqr3_ISW_constructT_check.o flops_geqr3_ISW_constructT.o
	$(LOADER) -o $@ main_flops_geqr3_ISW.o flops_geqr3_ISW_check.o flops_gemm.o flops_trmm.o flops_geqr3_check.o flops_lapack_larfg.o flops_geqr3_ISW.o flops_geqr3_ISW_save.o flops_geqr3.o flops_larft.o flops_lapack_geqr2.o flops_geqr3_bef_constructT_check.o flops_geqr3_bef_useT.o flops_geqr3_bef_useT_check.o flops_geqr3_wob_check.o flops_geqr3_ISW_constructT_check.o flops_geqr3_ISW_constructT.o

xmain_flops_geqr3.exe: main_flops_geqr3.o flops_geqr3_check.o flops_lapack_larfg.o flops_lapack_geqr2.o flops_geqr3_wob_check.o flops_geqr3_bef_constructT_check.o flops_geqr3_bef_useT_check.o  flops_geqr3_bef_constructT.o flops_geqr3_bef_useT.o flops_larft.o flops_geqr3.o
	$(LOADER) -o $@ main_flops_geqr3.o flops_geqr3_check.o flops_lapack_larfg.o flops_lapack_geqr2.o flops_geqr3_wob_check.o flops_geqr3_bef_constructT_check.o  flops_geqr3_bef_useT_check.o  flops_geqr3_bef_constructT.o  flops_geqr3_bef_useT.o  flops_larft.o flops_geqr3.o

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

xmain_flops_lapack_geqrf.exe: main_flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_geqrf.o flops_lapack_geqrf_check.o flops_lapack_larf.o flops_lapack_geqr2.o flops_larft.o flops_lapack_larfg.o flops_lapack_geqr2_check.o flops_lapack_geqrf_bef.o flops_lapack_geqrf_bef_check.o
	$(LOADER) -o $@  main_flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_geqrf.o flops_lapack_geqrf_check.o flops_lapack_larf.o flops_lapack_geqr2.o flops_larft.o flops_lapack_larfg.o flops_lapack_geqr2_check.o flops_lapack_geqrf_bef.o flops_lapack_geqrf_bef_check.o

xmain_flops_lapack_geqr2.exe: main_flops_lapack_geqr2.o flops_lapack_geqr2.o flops_lapack_geqr2_check.o flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_larfg.o flops_larft.o flops_lapack_larf.o flops_lapack_geqrf_check.o
	$(LOADER) -o $@  main_flops_lapack_geqr2.o flops_lapack_geqr2.o flops_lapack_geqr2_check.o flops_lapack_geqrf.o flops_lapack_larfb.o flops_lapack_larfg.o flops_larft.o flops_lapack_larf.o flops_lapack_geqrf_check.o

xmain_flops_lapack_orgqr.exe: main_flops_lapack_orgqr.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o flops_lapack_larfb.o flops_lapack_org2r.o flops_lapack_orgqr_bef_check.o flops_lapack_orgqr_bef.o
	$(LOADER) -o $@  main_flops_lapack_orgqr.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o  flops_lapack_larfb.o flops_lapack_org2r.o flops_lapack_orgqr_bef_check.o flops_lapack_orgqr_bef.o

xmain_flops_lapack_org2r.exe: main_flops_lapack_org2r.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o flops_lapack_larfb.o flops_lapack_org2r.o flops_lapack_org2r_Q1.o flops_lapack_org2r_Q1_check.o flops_lapack_org2r_Q2.o flops_lapack_org2r_Q2_check.o
	$(LOADER) -o $@  main_flops_lapack_org2r.o flops_lapack_orgqr_check.o flops_larft.o flops_lapack_org2r_check.o  flops_lapack_larfb.o flops_lapack_org2r.o flops_lapack_org2r_Q1.o flops_lapack_org2r_Q1_check.o flops_lapack_org2r_Q2.o flops_lapack_org2r_Q2_check.o

xmain_flops_legacy_lapack_org2r.exe: main_flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r_check.o flops_legacy_lapack_larf.o flops_lapack_larfb.o
	$(LOADER) -o $@  main_flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r.o flops_legacy_lapack_org2r_check.o flops_legacy_lapack_larf.o flops_lapack_larfb.o

xmain_flops_legacy_lapack_geqr2.exe: main_flops_legacy_lapack_geqr2.o flops_legacy_lapack_geqr2.o flops_legacy_lapack_geqr2_check.o flops_legacy_lapack_larf.o flops_legacy_lapack_larfg.o flops_lapack_larfb.o
	$(LOADER) -o $@  main_flops_legacy_lapack_geqr2.o flops_legacy_lapack_geqr2.o flops_legacy_lapack_geqr2_check.o flops_legacy_lapack_larf.o flops_legacy_lapack_larfg.o flops_lapack_larfb.o

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe *a
