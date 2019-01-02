
include make.inc

#all:  xmain_demo_02a.exe xmain_demo_03a.exe xmain_middle_02a.exe  xmain_unittest__lila_dgeqrf_w03_l
#all:  xmain_unittest__lila_dgeqrf_w03_l
all:  xmain_middle_02a.exe

lila_dgeqrf_w02a.o: lila_dgeqrf_w02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_w02a.c

lila_dgeqr2_w02a.o: lila_dgeqr2_w02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w02a.c

lila_wsq_dgeqr2_w02a.o: lila_wsq_dgeqr2_w02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_wsq_dgeqr2_w02a.c

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

lila_wsq_dormqrf_z02.o: lila_wsq_dormqrf_z02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_wsq_dormqrf_z02.c

lila_wsq_dormqrbz_w02.o: lila_wsq_dormqrbz_w02.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_wsq_dormqrbz_w02.c

lila_dgeqrf_w03a.o: lila_dgeqrf_w03a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqrf_w03a.c

lila_dgeqr2_w03a.o: lila_dgeqr2_w03a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dgeqr2_w03a.c

lila_dormqrbz_w03.o: lila_dormqrbz_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrbz_w03.c

lila_dormqrf_z03.o: lila_dormqrf_z03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_z03.c

lila_dormqrf_w03.o: lila_dormqrf_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dormqrf_w03.c

lila_dlarft_connect_w03.o: lila_dlarft_connect_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_connect_w03.c

lila_dlarft_w03.o: lila_dlarft_w03.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dlarft_w03.c

main_demo_02a.o: main_demo_02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_demo_02a.c

main_demo_03a.o: main_demo_03a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_demo_03a.c

main_unittest__lila_dgeqrf_w03_l.o: main_unittest__lila_dgeqrf_w03_l.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_unittest__lila_dgeqrf_w03_l.c

main_middle_02a.o: main_middle_02a.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_middle_02a.c

xmain_demo_02a.exe: main_demo_02a.o lila_dormqrf_w02.o lila_dgeqr2_w02a.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dormqrf_z00.o lila_dormqrf_z02.o lila_dlarft_connect_w02.o lila_wsq_dormqrf_z02.o lila_wsq_dgeqr2_w02a.o lila_wsq_dormqrbz_w02.o
	$(LOADER) -o $@ main_demo_02a.o lila_dormqrf_w02.o lila_dgeqr2_w02a.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dlarft_connect_w02.o lila_dormqrf_z00.o lila_dormqrf_z02.o lila_wsq_dormqrf_z02.o lila_wsq_dgeqr2_w02a.o lila_wsq_dormqrbz_w02.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_demo_03a.exe: main_demo_03a.o lila_dgeqrf_w03a.o lila_dormqrf_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o lila_dormqrf_z03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dormqrf_z00.o
	$(LOADER) -o $@ main_demo_03a.o lila_dgeqrf_w03a.o lila_dormqrf_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o lila_dormqrf_z03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dormqrf_z00.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_unittest__lila_dgeqrf_w03_l: main_unittest__lila_dgeqrf_w03_l.o lila_dgeqrf_w03a.o lila_dormqrf_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o lila_dormqrf_z03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dormqrf_z00.o
	$(LOADER) -o $@ main_unittest__lila_dgeqrf_w03_l.o lila_dgeqrf_w03a.o lila_dormqrf_w03.o lila_dgeqr2_w03a.o lila_dormqrbz_w03.o lila_dormqrf_z03.o lila_dlarft_connect_w03.o lila_dlarft_w03.o lila_dormqrf_z00.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

xmain_middle_02a.exe: main_middle_02a.o lila_dormqrf_w02.o lila_dgeqrf_w02a.o lila_dgeqr2_w02a.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dormqrf_z00.o lila_dormqrf_z02.o lila_dlarft_connect_w02.o lila_wsq_dormqrf_z02.o lila_wsq_dgeqr2_w02a.o lila_wsq_dormqrbz_w02.o
	$(LOADER) -o $@ main_middle_02a.o lila_dormqrf_w02.o lila_dgeqrf_w02a.o lila_dgeqr2_w02a.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dlarft_connect_w02.o lila_dormqrf_z00.o lila_dormqrf_z02.o lila_wsq_dormqrf_z02.o lila_wsq_dgeqr2_w02a.o lila_wsq_dormqrbz_w02.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS)

clean:
	rm -f lila_dlarft_connect_w02.o lila_dormqrbz_w02.o lila_dormqrf_w00.o lila_dormqrf_w02.o lila_dgeqr2_w02a.o lila_dormqrf_z00.o lila_dormqrf_z02.o main_demo_02a.o xmain_demo_02a.exe lila_wsq_dormqrf_z02.o lila_wsq_dgeqr2_w02a.o lila_wsq_dormqrbz_w02.o main_middle_02a.o lila_dgeqrf_w02a.o xmain_middle_02a.exe xmain_unittest__lila_dgeqrf_w03_l xmain_demo_03a.exe lila_dgeqr2_w03a.o lila_dlarft_connect_w03.o lila_dormqrbz_w03.o lila_dormqrf_z03.o main_unittest__lila_dgeqrf_w03_l.o lila_dgeqrf_w03a.o lila_dlarft_w03.o lila_dormqrf_w03.o main_demo_03a.o

