
include ../make.inc

all:  xmain_LAPACK.exe  

xN2T.o: xN2T.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) xN2T.c

xV2N.o: xV2N.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) xV2N.c

nV2T.o: nV2T.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) nV2T.c

main_LAPACK.o: main_LAPACK.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_LAPACK.c

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

xmain_LAPACK.exe: main_LAPACK.o liblila.a nV2T.o xN2T.o xV2N.o
	$(LOADER) -o $@ main_LAPACK.o ../src/liblila.a nV2T.o xN2T.o xV2N.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

clean:
	rm -f  xmain_LAPACK.exe main_LAPACK.o nV2T.o xN2T.o xV2N.o
 
