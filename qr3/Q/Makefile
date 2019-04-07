
include ../make.inc

all: xmain_dorgqr.exe xmain_flops_V2N.exe xmain_flops_org2r.exe

xmain_dorgqr.exe: main_dorgqr.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_larfb.o flops_N2T.o flops_VT2Q.o 
	$(LOADER) -o $@ main_dorgqr.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dorgqr.o dV2Q.o dV2N.o dN2T.o dVT2Q.o ULTinU.o mLUinA.o ApUBTinA.o flops_org2r.o flops_gemm.o flops_trmm.o flops_syrk.o flops_V2N.o flops_larft.o flops_larfb.o flops_N2T.o flops_VT2Q.o  $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_flops_V2N.exe: main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o
	$(LOADER) -o $@ main_flops_V2N.o flops_trmm.o flops_syrk.o flops_V2N.o

xmain_flops_org2r.exe: main_flops_org2r.o flops_org2r.o
	$(LOADER) -o $@ main_flops_org2r.o flops_org2r.o

.c.o:
	$(CC) $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -c -o $@ $<

clean:
	rm -f *o *exe
