
include make.inc

all:  xmain_dgmres.exe 

matvec_Mr.o: matvec_Mr.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) matvec_Mr.c

matvec_Ml.o: matvec_Ml.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) matvec_Ml.c

matvec_A.o: matvec_A.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) matvec_A.c

dgmres.o: dgmres.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) dgmres.c

main_dgmres.o: main_dgmres.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) -I$(INCLUDELAPACKE) main_dgmres.c

xmain_dgmres.exe: main_dgmres.o matvec_A.o matvec_Ml.o matvec_Mr.o dgmres.o
	$(LOADER) -o $@ main_dgmres.o matvec_A.o matvec_Ml.o matvec_Mr.o dgmres.o $(LIBLAPACKE) $(LIBLAPACK) $(LIBCBLAS) $(LIBBLAS) $(LIBS) 

clean:
	rm -f main_dgmres.o matvec_A.o matvec_Ml.o matvec_Mr.o dgmres.o

