
include make.inc

all:  xmain.exe

lila_dge_qr_larft_connect_w02.o: lila_dge_qr_larft_connect_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_larft_connect_w02.c

lila_dge_qr_larft_connect_w03.o: lila_dge_qr_larft_connect_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_larft_connect_w03.c

lila_dge_qr_ormqrbz_w02.o: lila_dge_qr_ormqrbz_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrbz_w02.c

lila_dge_qr_ormqrbz_w00.o: lila_dge_qr_ormqrbz_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrbz_w00.c

lila_dge_qr_ormqrbz_vL0.o: lila_dge_qr_ormqrbz_vL0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrbz_vL0.c

lila_dge_qr_ormqrf_vL0.o: lila_dge_qr_ormqrf_vL0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrf_vL0.c

lila_dge_qr_ormqrf_w00.o: lila_dge_qr_ormqrf_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrf_w00.c

lila_dge_qr_ormqrf_w02.o: lila_dge_qr_ormqrf_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrf_w02.c

lila_dge_qr_ormqrf_w03.o: lila_dge_qr_ormqrf_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_ormqrf_w03.c

lila_dge_qr_wq_vr0.o: lila_dge_qr_wq_vr0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_vr0.c

lila_dge_qr_wq_w02.o: lila_dge_qr_wq_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_w02.c

lila_dge_qr_wq_vL0.o: lila_dge_qr_wq_vL0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_vL0.c

lila_dge_qr_wq_manylevels_vr0.o: lila_dge_qr_wq_manylevels_vr0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_manylevels_vr0.c

lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0.o: lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0.c

lila_dge_qr_wq_levelx_w00.o: lila_dge_qr_wq_levelx_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_levelx_w00.c

lila_dge_qr_wq_levelx_w02.o: lila_dge_qr_wq_levelx_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_levelx_w02.c

lila_dge_qr_wq_levelx_w03.o: lila_dge_qr_wq_levelx_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_levelx_w03.c

lila_dlarft_w03.o: lila_dlarft_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_w03.c

dgeqr3.o: dgeqr3.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) dgeqr3.c

main.o: main.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main.c

xmain.exe: main.o dgeqr3.o lila_dge_qr_wq_vL0.o lila_dge_qr_wq_vr0.o lila_dge_qr_ormqrf_vL0.o lila_dge_qr_ormqrbz_vL0.o lila_dge_qr_wq_manylevels_vr0.o lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0.o lila_dge_qr_wq_levelx_w00.o lila_dge_qr_wq_levelx_w02.o lila_dge_qr_ormqrf_w02.o lila_dge_qr_ormqrbz_w02.o lila_dge_qr_wq_w02.o lila_dge_qr_larft_connect_w02.o lila_dge_qr_ormqrf_w00.o lila_dge_qr_ormqrbz_w00.o lila_dge_qr_wq_levelx_w03.o lila_dge_qr_ormqrf_w03.o lila_dge_qr_larft_connect_w03.o lila_dlarft_w03.o
	$(LOADER) -o $@ main.o dgeqr3.o lila_dge_qr_wq_vL0.o lila_dge_qr_wq_vr0.o lila_dge_qr_ormqrf_vL0.o lila_dge_qr_ormqrbz_vL0.o lila_dge_qr_wq_manylevels_vr0.o lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0.o lila_dge_qr_wq_levelx_w00.o lila_dge_qr_wq_levelx_w02.o lila_dge_qr_ormqrf_w02.o lila_dge_qr_ormqrbz_w02.o lila_dge_qr_wq_w02.o lila_dge_qr_larft_connect_w02.o lila_dge_qr_ormqrf_w00.o lila_dge_qr_ormqrbz_w00.o lila_dge_qr_wq_levelx_w03.o lila_dge_qr_ormqrf_w03.o lila_dge_qr_larft_connect_w03.o lila_dlarft_w03.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS) 

clean:
	rm -f xmain.exe main.o dgeqr3.o lila_dge_qr_wq_vL0.o lila_dge_qr_wq_vr0.o lila_dge_qr_ormqrf_vL0.o lila_dge_qr_ormqrbz_vL0.o lila_dge_qr_wq_manylevels_vr0.o lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0.o lila_dge_qr_wq_levelx_w00.o lila_dge_qr_wq_levelx_w02.o lila_dge_qr_ormqrf_w02.o lila_dge_qr_ormqrbz_w02.o lila_dge_qr_wq_w02.o lila_dge_qr_larft_connect_w02.o lila_dge_qr_ormqrf_w00.o lila_dge_qr_ormqrbz_w00.o lila_dge_qr_wq_levelx_w03.o lila_dge_qr_ormqrf_w03.o lila_dge_qr_larft_connect_w03.o lila_dlarft_w03.o

