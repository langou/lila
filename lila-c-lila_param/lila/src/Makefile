
include ../make.inc

SRC = lila_dgeqrf.o dgeqr3.o lila_dgeqrf_recursive.o lila_dgeqrf_levelx.o lila_dormqrf.o lila_dlarft_connect.o lila_dormqrbz.o lila_dlarft.o \
	lila_dormqrf_z03.o \
	lila_dorgh2_3.o lila_ormhr2_3.o lila_ormhr2_w03_hr.o \
	lila_dgeqrf_mt.o lila_dgeqrf_mt_l.o lila_dgeqrf_mt_hr.o lila_dgeqrf_mt_hh.o \
	lila_dgeqr2.o lila_dgeqr2_l.o lila_dgeqr2_3.o lila_dgeqr2_hr.o \
	lila_dT2tau.o lila_dV2tau_w03.o \
	lila_query_dormqrbz_w03.o lila_query_dgeqrf_w03_mt_hr.o lila_query_dormqrf_w03.o lila_query_dgeqrf_w03_mt.o \
	lila_query_dgeqrf_w03_mt_l.o lila_query_dgeqr2_w03_hr.o lila_query_dgeqr2_w03_3.o lila_query_dgeqr2_w03_l.o \
	lila_query_ormhr2_w03_hr.o lila_query_dorgh2_3.o lila_query_dgeqrf_w03_levelx.o lila_query_dgeqrf_w03_recursive.o \
	lila_query_dgeqrf_w03_mt_hh.o \
	lila_test_qr_repres_1.o lila_test_qq_orth_1.o lila_test_r_repres_2.o lila_test_hh_repres.o lila_test_vt_repres.o

CC = gcc
OPTS =

all: $(SRC)
	ar -cr liblila.a $^
	ranlib liblila.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) -I$(INCLUDELILA) -c -o $@ $<

clean:
	rm -f *o *a

