
include ../make.inc

FLOPS = \
	flops_geqr3.o         \
	flops_lapack_geqr2.o  \
	flops_lapack_org2r.o  \
	flops_larft.o

FLOPS_RESEARCH = \
	flops_geqr3.o                      \
	flops_geqr3_check.o                \
	flops_geqr3_wob_check.o            \
	flops_geqr3_bef_constructT.o       \
	flops_geqr3_bef_constructT_check.o \
	flops_geqr3_bef_useT.o             \
	flops_geqr3_bef_useT_check.o       \
	flops_geqr3_ISW_constructT.o       \
	flops_geqr3_ISW_constructT_check.o \
	flops_lapack_larfg.o               \
	flops_lapack_larfb.o               \
	flops_lapack_larfb_check.o         \
	flops_lapack_geqr2.o               \
	flops_lapack_geqr2_check.o         \
	flops_lapack_geqrf.o               \
	flops_lapack_geqrf_check.o         \
	flops_lapack_larf.o                \
	flops_lapack_geqrf_bef.o           \
	flops_lapack_geqrf_bef_check.o     \
	flops_lapack_orgqr_check.o         \
	flops_lapack_org2r.o               \
	flops_lapack_org2r_check.o         \
	flops_lapack_orgqr_bef.o           \
	flops_lapack_orgqr_bef_check.o 	   \
	flops_lapack_org2r_Q1.o            \
	flops_lapack_org2r_Q1_check.o      \
	flops_lapack_org2r_Q2.o            \
	flops_lapack_org2r_Q2_check.o      \
	flops_legacy_lapack_geqr2.o        \
	flops_legacy_lapack_geqr2_check.o  \
	flops_legacy_lapack_org2r.o        \
	flops_legacy_lapack_org2r_check.o  \
	flops_legacy_lapack_larfg.o        \
	flops_legacy_lapack_larf.o         \
	flops_larft.o                      \
	flops_gemm.o                       \
	flops_trmm.o                       \
	flops_syrk.o                       \
	flops_N2T.o                        \
	flops_N2T_check.o                  \
	flops_ULTinU.o                     \
	flops_mLUinA.o                     \
 	flops_VT2Q.o                       \
	flops_VT2Q_check.o                 \
	flops_V2N.o                        \
	flops_V2N_check.o                  \
	flops_mLUinA_check.o               \
	flops_ApUBTinA.o                   \
	flops_ApUBTinA_check.o             \
	flops_ULTinU.o                     \
	flops_ULTinU_check.o               \
	flops_geqr3_UT_check.o             \
	flops_geqr3_UT_save_check.o        \
	flops_geqr3_UT.o                   \
	flops_geqr3_ISW_check.o            \
	flops_geqr3_ISW_UT.o               \
	flops_geqr3_ISW_UT_check.o         \
	flops_geqr3_ISW.o

all: xmain_flops_lapack_orgqr.exe xmain_flops_lapack_org2r.exe xmain_flops_lapack_geqr2.exe xmain_flops_lapack_geqrf.exe xmain_flops_geqr3.exe xmain_flops_geqr3_ISW.exe xmain_flops_N2T.exe xmain_flops_V2N.exe xmain_flops_VT2Q.exe xmain_flops_mLUinA.exe xmain_flops_ApUBTinA.exe xmain_flops_ULTinU.exe xmain_flops_lapack_larfb.exe xmain_flops_lapack_geqrf.exe xmain_flops_lapack_geqr2.exe xmain_flops_lapack_orgqr.exe xmain_flops_lapack_org2r.exe xmain_flops_legacy_lapack_org2r.exe xmain_flops_legacy_lapack_geqr2.exe xmain_flops_geqr3_UT.exe xmain_flops_geqr3_ISW_UT.exe

#all: xmain_flops_geqr3.exe 
#all: xmain_flops_geqr3_ISW.exe

libflops.a: $(FLOPS) $(FLOPS_RESEARCH)
	ar -cr libflops.a $^
	ranlib libflops.a

xmain_flops_geqr3_ISW_UT.exe: main_flops_geqr3_ISW_UT.o libflops.a
	$(LOADER) -o $@ main_flops_geqr3_ISW_UT.o libflops.a

xmain_flops_geqr3_UT.exe: main_flops_geqr3_UT.o libflops.a
	$(LOADER) -o $@ main_flops_geqr3_UT.o libflops.a

xmain_flops_geqr3_ISW.exe: main_flops_geqr3_ISW.o libflops.a
	$(LOADER) -o $@ main_flops_geqr3_ISW.o libflops.a

xmain_flops_geqr3.exe: main_flops_geqr3.o libflops.a
	$(LOADER) -o $@ main_flops_geqr3.o libflops.a

xmain_flops_N2T.exe: main_flops_N2T.o libflops.a 
	$(LOADER) -o $@ main_flops_N2T.o libflops.a

xmain_flops_V2N.exe: main_flops_V2N.o libflops.a 
	$(LOADER) -o $@ main_flops_V2N.o libflops.a 

xmain_flops_VT2Q.exe: main_flops_VT2Q.o libflops.a 
	$(LOADER) -o $@ main_flops_VT2Q.o libflops.a 

xmain_flops_mLUinA.exe: main_flops_mLUinA.o libflops.a 
	$(LOADER) -o $@ main_flops_mLUinA.o libflops.a

xmain_flops_ApUBTinA.exe: main_flops_ApUBTinA.o libflops.a
	$(LOADER) -o $@  main_flops_ApUBTinA.o libflops.a

xmain_flops_ULTinU.exe: main_flops_ULTinU.o libflops.a
	$(LOADER) -o $@  main_flops_ULTinU.o libflops.a

xmain_flops_lapack_larfb.exe: main_flops_lapack_larfb.o libflops.a 
	$(LOADER) -o $@  main_flops_lapack_larfb.o libflops.a

xmain_flops_lapack_geqrf.exe: main_flops_lapack_geqrf.o libflops.a
	$(LOADER) -o $@  main_flops_lapack_geqrf.o libflops.a

xmain_flops_lapack_geqr2.exe: main_flops_lapack_geqr2.o libflops.a
	$(LOADER) -o $@  main_flops_lapack_geqr2.o libflops.a

xmain_flops_lapack_orgqr.exe: main_flops_lapack_orgqr.o libflops.a
	$(LOADER) -o $@  main_flops_lapack_orgqr.o libflops.a

xmain_flops_lapack_org2r.exe: main_flops_lapack_org2r.o libflops.a
	$(LOADER) -o $@  main_flops_lapack_org2r.o libflops.a

xmain_flops_legacy_lapack_org2r.exe: main_flops_legacy_lapack_org2r.o libflops.a
	$(LOADER) -o $@  main_flops_legacy_lapack_org2r.o libflops.a

xmain_flops_legacy_lapack_geqr2.exe: main_flops_legacy_lapack_geqr2.o libflops.a
	$(LOADER) -o $@  main_flops_legacy_lapack_geqr2.o libflops.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe *a
