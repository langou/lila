
include ../make.inc

include ../make.inc

all: xmain_VT2Q.exe 
all: xmain_dorgqr.exe 

xmain_VT2Q.exe: main_VT2Q.o VT2Q.o ULTinU.o mLUinA.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o ApUBTinA.o
	$(LOADER) -o $@ main_VT2Q.o VT2Q.o ULTinU.o mLUinA.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o ApUBTinA.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_dorgqr.exe: main_dorgqr.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dorgqr.o xV2N.o xN2T.o
	$(LOADER) -o $@ main_dorgqr.o qr3_test_qq_orth_1.o qr3_test_qr_repres_1.o dorgqr.o xV2N.o xN2T.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) -c -o $@ $<

clean:
	rm -f *o *exe
