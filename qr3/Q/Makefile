
include ../make.inc

all:  xmain_VT2Q.exe  

VT2Q.o: VT2Q.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) VT2Q.c

ULTinU.o: ULTinU.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) ULTinU.c

LUinQ.o: LUinQ.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) LUinQ.c

main_VT2Q.o: main_VT2Q.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) main_VT2Q.c

liblila.a: ../src/dgeqr3.c \
	../src/lila_dgeqrf_mt.c \
	../src/lila_dormqrbz.c \
	../src/lila_query_dgeqr2_w03_l.c \
	../src/lila_query_dormqrbz_w03.c \
	../src/lila_dT2tau.c \
	../src/lila_dgeqrf_mt_hh.c \
	../src/lila_dormqrf.c \
	../src/lila_query_dgeqrf_w03_levelx.c \
	../src/lila_query_dormqrf_w03.c \
	../src/lila_dV2tau_w03.c \
	../src/lila_dgeqrf_mt_hr.c \
	../src/lila_query_dgeqrf_w03_mt.c \
	../src/lila_query_ormhr2_w03_hr.c \
	../src/lila_dgeqr2.c \
	../src/lila_dgeqrf_mt_l.c \
	../src/lila_dormqrf_z03.c \
	../src/lila_query_dgeqrf_w03_mt_hh.c \
	../src/lila_test_hh_repres.c \
	../src/lila_test_vt_repres.c \
	../src/lila_dgeqr2_3.c \
	../src/lila_dgeqrf_recursive.c \
	../src/lila_ormhr2_3.c \
	../src/lila_query_dgeqrf_w03_mt_hr.c \
	../src/lila_test_qq_orth_1.c \
	../src/lila_dgeqr2_hr.c \
	../src/lila_dlarft.c \
	../src/lila_ormhr2_w03_hr.c \
	../src/lila_query_dgeqrf_w03_mt_l.c \
	../src/lila_test_qr_repres_1.c \
	../src/lila_dgeqr2_l.c \
	../src/lila_dlarft_connect.c \
	../src/lila_query_dgeqr2_w03_3.c \
	../src/lila_query_dgeqrf_w03_recursive.c \
	../src/lila_test_r_repres_2.c \
	../src/lila_dgeqrf.c \
	../src/lila_dorgh2_3.c \
	../src/lila_query_dgeqr2_w03_hr.c \
	../src/lila_query_dorgh2_3.c \
	../src/lila_dgeqrf_levelx.c 
	make -C ../src/

xmain_VT2Q.exe: main_VT2Q.o liblila.a VT2Q.o ULTinU.o LUinQ.o
	$(LOADER) -o $@ main_VT2Q.o ../src/liblila.a VT2Q.o ULTinU.o LUinQ.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

clean:
	rm -f  xmain_VT2Q.exe main_VT2Q.o VT2Q.o ULTinU.o LUinQ.o
 
