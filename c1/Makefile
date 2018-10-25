
include make.inc

all:  xmain.exe

lila_dlarft_connect_w02.o: lila_dlarft_connect_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_connect_w02.c

lila_dlarft_connect_w03.o: lila_dlarft_connect_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_connect_w03.c

lila_dormqrbz_w03.o: lila_dormqrbz_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrbz_w03.c

lila_dormqrbz_w02.o: lila_dormqrbz_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrbz_w02.c

lila_dormqrbz_w00.o: lila_dormqrbz_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrbz_w00.c

lila_dormqrf_w00.o: lila_dormqrf_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_w00.c

lila_dormqrf_w02.o: lila_dormqrf_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_w02.c

lila_dormqrf_w03.o: lila_dormqrf_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_w03.c

lila_dgeqr2_w02a.o: lila_dgeqr2_w02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w02a.c

lila_dgeqr2_w02b.o: lila_dgeqr2_w02b.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w02b.c

lila_dgeqr2_w03b.o: lila_dgeqr2_w03b.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w03b.c

lila_dgeqr2_w03a.o: lila_dgeqr2_w03a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w03a.c

lila_dgeqrf_levelx_w00.o: lila_dgeqrf_levelx_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_levelx_w00.c

lila_dgeqrf_levelx_w02.o: lila_dgeqrf_levelx_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_levelx_w02.c

lila_dgeqrf_levelx_w03.o: lila_dgeqrf_levelx_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_levelx_w03.c

lila_dgeqrf_recursive_w03.o: lila_dgeqrf_recursive_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_recursive_w03.c

lila_dgeqrf_recursive_w02.o: lila_dgeqrf_recursive_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_recursive_w02.c

lila_dlarft_w03.o: lila_dlarft_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_w03.c

lila_dgeqr1_w03a.o: lila_dgeqr1_w03a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr1_w03a.c

lila_dgeqr1_w03b.o: lila_dgeqr1_w03b.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr1_w03b.c

dgeqr3.o: dgeqr3.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) dgeqr3.c

main.o: main.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main.c

xmain.exe: main.o dgeqr3.o lila_dgeqrf_levelx_w00.o lila_dgeqrf_levelx_w02.o lila_dormqrf_w02.o lila_dormqrbz_w02.o lila_dgeqr2_w02a.o lila_dlarft_connect_w02.o lila_dormqrf_w00.o lila_dormqrbz_w00.o lila_dgeqrf_levelx_w03.o lila_dormqrf_w03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o  lila_dgeqr2_w02b.o lila_dgeqrf_recursive_w02.o lila_dgeqrf_recursive_w03.o lila_dgeqr1_w03a.o lila_dgeqr1_w03b.o lila_dgeqr2_w03b.o
	$(LOADER) -o $@ main.o dgeqr3.o lila_dgeqrf_levelx_w00.o lila_dgeqrf_levelx_w02.o lila_dormqrf_w02.o lila_dormqrbz_w02.o lila_dgeqr2_w02a.o lila_dlarft_connect_w02.o lila_dormqrf_w00.o lila_dormqrbz_w00.o lila_dgeqrf_levelx_w03.o lila_dormqrf_w03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o lila_dgeqr2_w02b.o lila_dgeqrf_recursive_w02.o lila_dgeqrf_recursive_w03.o lila_dgeqr1_w03a.o lila_dgeqr1_w03b.o lila_dgeqr2_w03b.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS) 

clean:
	rm -f xmain.exe main.o dgeqr3.o lila_dgeqrf_levelx_w00.o lila_dgeqrf_levelx_w02.o lila_dormqrf_w02.o lila_dormqrbz_w02.o lila_dgeqr2_w02a.o lila_dlarft_connect_w02.o lila_dormqrf_w00.o lila_dormqrbz_w00.o lila_dgeqrf_levelx_w03.o lila_dormqrf_w03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o lila_dgeqr2_w02b.o lila_dgeqrf_recursive_w02.o lila_dgeqrf_recursive_w03.o lila_dgeqr1_w03a.o lila_dgeqr1_w03b.o lila_dgeqr2_w03b.o
