
include ../make.inc

#all: xmain_flops_lapack_orgqr.exe xmain_flops_lapack_org2r.exe xmain_flops_lapack_geqr2.exe xmain_flops_lapack_geqrf.exe xmain_flops_geqr3.exe xmain_flops_geqr3_ISW.exe xmain_flops_N2T.exe xmain_flops_V2N.exe xmain_flops_VT2Q.exe xmain_flops_mLUinA.exe xmain_flops_ApUBTinA.exe xmain_flops_ULTinU.exe xmain_flops_lapack_larfb.exe xmain_flops_lapack_geqrf.exe xmain_flops_lapack_geqr2.exe xmain_flops_lapack_orgqr.exe xmain_flops_lapack_org2r.exe xmain_flops_legacy_lapack_org2r.exe xmain_flops_legacy_lapack_geqr2.exe xmain_flops_geqr3_UT.exe xmain_flops_geqr3_ISW_UT.exe xmain_flops_qr2_dorgqr.exe xmain_flops_larft3.exe

#all: xmain_flops_geqr3.exe 
#all: xmain_flops_geqr3_ISW.exe
#all: xmain_flops_lapack_org2r.exe
#all: xmain_flops_qr2_dorgqr.exe
all: xmain_flops_geqr3_noR.exe
#all: xmain_flops_larft3.exe // Same as flops_larft

libflops.a: $(FLOPS) $(FLOPS_RESEARCH)
	ar -cr libflops.a $^
	ranlib libflops.a

xmain_flops_larft3.exe: main_flops_larft3.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_larft3.o ../flops/libflops.a $(LIB)

xmain_flops_qr2_dorgqr.exe: main_flops_qr2_dorgqr.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_qr2_dorgqr.o ../flops/libflops.a $(LIB)

xmain_flops_geqr3_ISW_UT.exe: main_flops_geqr3_ISW_UT.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_geqr3_ISW_UT.o ../flops/libflops.a $(LIB)

xmain_flops_geqr3_UT.exe: main_flops_geqr3_UT.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_geqr3_UT.o ../flops/libflops.a $(LIB)

xmain_flops_geqr3_ISW.exe: main_flops_geqr3_ISW.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_geqr3_ISW.o ../flops/libflops.a $(LIB)

xmain_flops_geqr3_noR.exe: main_flops_geqr3_noR.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_geqr3_noR.o ../flops/libflops.a $(LIB)

xmain_flops_geqr3.exe: main_flops_geqr3.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_geqr3.o ../flops/libflops.a $(LIB)

xmain_flops_N2T.exe: main_flops_N2T.o  
	make -C ../flops/
	$(LOADER) -o $@ main_flops_N2T.o ../flops/libflops.a $(LIB)

xmain_flops_V2N.exe: main_flops_V2N.o  
	make -C ../flops/
	$(LOADER) -o $@ main_flops_V2N.o ../flops/libflops.a $(LIB)

xmain_flops_VT2Q.exe: main_flops_VT2Q.o  
	make -C ../flops/
	$(LOADER) -o $@ main_flops_VT2Q.o ../flops/libflops.a $(LIB)

xmain_flops_mLUinA.exe: main_flops_mLUinA.o  
	make -C ../flops/
	$(LOADER) -o $@ main_flops_mLUinA.o ../flops/libflops.a $(LIB)

xmain_flops_ApUBTinA.exe: main_flops_ApUBTinA.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_ApUBTinA.o ../flops/libflops.a $(LIB)

xmain_flops_ULTinU.exe: main_flops_ULTinU.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_ULTinU.o ../flops/libflops.a $(LIB)

xmain_flops_lapack_larfb.exe: main_flops_lapack_larfb.o  
	make -C ../flops/
	$(LOADER) -o $@ main_flops_lapack_larfb.o ../flops/libflops.a $(LIB)

xmain_flops_lapack_geqrf.exe: main_flops_lapack_geqrf.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_lapack_geqrf.o ../flops/libflops.a $(LIB)

xmain_flops_lapack_geqr2.exe: main_flops_lapack_geqr2.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_lapack_geqr2.o ../flops/libflops.a $(LIB)

xmain_flops_lapack_orgqr.exe: main_flops_lapack_orgqr.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_lapack_orgqr.o ../flops/libflops.a $(LIB)

xmain_flops_lapack_org2r.exe: main_flops_lapack_org2r.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_lapack_org2r.o ../flops/libflops.a $(LIB)

xmain_flops_legacy_lapack_org2r.exe: main_flops_legacy_lapack_org2r.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_legacy_lapack_org2r.o ../flops/libflops.a $(LIB)

xmain_flops_legacy_lapack_geqr2.exe: main_flops_legacy_lapack_geqr2.o 
	make -C ../flops/
	$(LOADER) -o $@ main_flops_legacy_lapack_geqr2.o ../flops/libflops.a $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe *a
