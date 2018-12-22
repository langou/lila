
include make.inc

all:  xmain_demo_02a.exe

lila_dgeqr2_w02a.o: lila_dgeqr2_w02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w02a.c

lila_dlarft_connect_w02.o: lila_dlarft_connect_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_connect_w02.c

lila_dormqrbz_w02.o: lila_dormqrbz_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrbz_w02.c

lila_dormqrf_w00.o: lila_dormqrf_w00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_w00.c

lila_dormqrf_w02.o: lila_dormqrf_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_w02.c

lila_dormqrf_z00.o: lila_dormqrf_z00.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_z00.c

lila_dormqrf_z02.o: lila_dormqrf_z02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_z02.c

main_demo_02a.o: main_demo_02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_demo_02a.c

xmain_demo_02a.exe: main_demo_02a.o lila_dormqrf_w02.o lila_dgeqr2_w02a.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dormqrf_z00.o lila_dormqrf_z02.o lila_dlarft_connect_w02.o
	$(LOADER) -o $@ main_demo_02a.o lila_dormqrf_w02.o lila_dgeqr2_w02a.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dlarft_connect_w02.o lila_dormqrf_z00.o lila_dormqrf_z02.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

clean:
	rm -f lila_dlarft_connect_w02.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dormqrf_w02.o lila_dgeqr2_w02a.o lila_dormqrf_z00.o lila_dormqrf_z02.o main_demo_02a.o xmain_demo_02a.exe 



