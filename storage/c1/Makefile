
include make.inc

all:  xmain.exe

lila_dge_qr_wq_vr0.o: lila_dge_qr_wq_vr0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_vr0.c

lila_dge_qr_wq_vL0.o: lila_dge_qr_wq_vL0.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) lila_dge_qr_wq_vL0.c

dgeqr3.o: dgeqr3.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) dgeqr3.c

main.o: main.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main.c

xmain.exe: main.o dgeqr3.o lila_dge_qr_wq_vL0.o lila_dge_qr_wq_vr0.o
	$(LOADER) -o $@ main.o dgeqr3.o lila_dge_qr_wq_vL0.o lila_dge_qr_wq_vr0.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS) 

clean:
	rm -f xmain.exe main.o dgeqr3.o lila_dge_qr_wq_vL0.o lila_dge_qr_wq_vr0.o

