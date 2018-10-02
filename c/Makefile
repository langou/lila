
include make.inc

all:  xmain.exe

dgeqr3.o: dgeqr3.c
	$(CC) -c $(OPTS)  -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE)  dgeqr3.c

main.o: main.c
	$(CC) -c $(OPTS)  -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE)  main.c

xmain.exe: main.o dgeqr3.o
	$(LOADER) -o $@ main.o dgeqr3.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS) 

clean:
	rm -f xmain.exe main.o dgeqr3.o

