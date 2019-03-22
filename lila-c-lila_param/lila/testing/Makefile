
include ../make.inc

all:  xmain_recursive.exe xmain_levelx.exe xmain_append_testing.exe

lila_query_dgeqrf_ker_recursive.o: lila_query_dgeqrf_ker_recursive.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_query_dgeqrf_ker_recursive.c

lila_query_dgeqrf_LAPACK_appendcols.o: lila_query_dgeqrf_LAPACK_appendcols.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_query_dgeqrf_LAPACK_appendcols.c

lila_dgeqrf_qr2_recursive.o: lila_dgeqrf_qr2_recursive.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_dgeqrf_qr2_recursive.c

lila_dgeqrf_qr2_recursive2.o: lila_dgeqrf_qr2_recursive2.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_dgeqrf_qr2_recursive2.c

lila_dgeqrf_qr2.o: lila_dgeqrf_qr2.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_dgeqrf_qr2.c

lila_dgeqrf_ker_recursive.o: lila_dgeqrf_ker_recursive.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_dgeqrf_ker_recursive.c

lila_dgeqrf_w03_appendcols.o: lila_dgeqrf_w03_appendcols.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_dgeqrf_w03_appendcols.c

lila_dgeqrf_LAPACK_appendcols.o: lila_dgeqrf_LAPACK_appendcols.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) lila_dgeqrf_LAPACK_appendcols.c

main_recursive.o: main_recursive.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) main_recursive.c

main_levelx.o: main_levelx.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) main_levelx.c

main_append_testing.o: main_append_testing.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) main_append_testing.c

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

xmain_recursive.exe: main_recursive.o liblila.a 
	$(LOADER) -o $@ main_recursive.o ../src/liblila.a $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_levelx.exe: main_levelx.o liblila.a 
	$(LOADER) -o $@ main_levelx.o ../src/liblila.a $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_append_testing.exe: main_append_testing.o liblila.a 
	$(LOADER) -o $@ main_append_testing.o ../src/liblila.a $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

clean:
	rm -f lila_dormqrf_w00.o lila_dormqrf_z00.o lila_dormqrf_z02.o main_demo_02a.o xmain_demo_02a.exe main_middle_02.o lila_dgeqrf_w02a.o xmain_middle_02.exe xmain_unittest__lila_dgeqrf_w03_l xmain_demo_03a.exe lila_dgeqrf_w03_mt_l.o lila_dlarft_connect_w03.o lila_dormqrbz_w03.o lila_dormqrf_z03.o main_unittest__lila_dgeqrf_w03_l.o lila_dgeqrf_w03a.o lila_dlarft_w03.o lila_dormqrf_w03.o main_demo_03a.o lila_dgeqrf_w03_level1.o lila_dgeqrf_w03_mt.o lila_dgeqr2_w03_l.o lila_dgeqrf_w03_l.o lila_dgeqr2_w03_3.o lila_dgeqrf_w02b.o lila_ormhr_w0b.o lila_dorgh2_3.o  lila_dgeqr2_w03_hr.o lila_ormhr2_3.o lila_ormhr2_w03_hr.o lila_dgeqrf_w03_mt_hr.o main_w03.o lila_dgeqrf_w03_levelx.o lila_dgeqrf_w03_recursive.o lila_query_dormqrbz_w03.o lila_query_dgeqrf_w03_mt_hr.o lila_query_dormqrf_w03.o lila_query_dgeqrf_w03_mt.o lila_query_dgeqrf_w03_mt_l.o lila_query_dgeqr2_w03_hr.o lila_query_dgeqr2_w03_3.o lila_query_dgeqr2_w03_l.o lila_query_ormhr2_w03_hr.o lila_query_dorgh2_3.o lila_query_dgeqrf_w03_recursive.o lila_query_dgeqrf_w03_levelx.o main_vq_LAPACK.o main_qr3.o xmain_vq_LAPACK.exe xmain_qr3.exe xmain_w03.exe main_v03.o xmain_v03.exe lila_dgeqrf_v03_levelx.o lila_dgeqrf_v03_mt.o lila_dgeqrf_v03_mt_l.o lila_dgeqr2_v03_l.o lila_dgeqr2_v03_3.o lila_dgeqrf_v03_recursive.o lila_dgeqrf_q03_3.o main_vq03.o xmain_vq03.exe lila_dgeqrf_q03_mt.o lila_dgeqrf_q03_mt_l.o lila_dgeqr2_q03_l.o lila_dgeqr2_q03_3.o lila_dgeqrf_q03_recursive.o lila_dgeqrf_q03_levelx.o lila_dT2tau_w03.o lila_dgeqrf_w03_mt_hh.o lila_dgeqr2_w03.o lila_query_dgeqrf_w03_mt_hh.o lila_main_test.o lila_dgeqrf_v03_mt_hh.o lila_dgeqrf_q03_mt_hh.o lila_dgeqr2_q03.o lila_dgeqr2_v03.o main_ker_LAPACK.o xmain_ker_LAPACK.exe xmain_qr2.exe main_qr2.o lila_dgeqrf_qr2_recursive.o lila_dgeqrf_qr2.o lila_dgeqrf_qr2_recursive2.o lila_dgeqrf_ker_recursive.o lila_query_dgeqrf_ker_recursive.o lila_dV2tau_w03.o main_w03_appendcols.o xmain_w03_appendcols.exe main_LAPACK_appendcols.o xmain_LAPACK_appendcols.exe lila_dgeqrf_w03_appendcols.o lila_dgeqrf_LAPACK_appendcols.o lila_query_dgeqrf_LAPACK_appendcols.o lila_test_qr_repres_1.o lila_dgeqrf_vrtq_levelx.o xmain_vrtq.exe lila_dgeqrf.o lila_dgeqrf_vrtq_recursive.o main_vrtq.o lila_dgeqrf_t03.o lila_dgeqr2_t03_l.o lila_dgeqrf_t03_l.o lila_test_qq_orth_1.o xmain.exe main.o lila_dgeqrf_recursive.o lila_dgeqrf_mt.o lila_dormqrf.o lila_dlarft_connect.o lila_dormqrbz.o lila_dgeqrf_mt_l.o lila_dgeqrf_mt_hr.o lila_dgeqrf_mt_hh.o lila_dgeqr2.o lila_dgeqr2_l.o lila_dgeqr2_3.o lila_dgeqr2_hr.o lila_dlarft.o lila_test_hh_repres.o xmain_recursive.o xmain_recursive.exe main_levelx.o xmain_levelx.exe main_append_testing.o xmain_append_testing.exe
 
