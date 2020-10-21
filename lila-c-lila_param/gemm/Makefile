
include make.inc

all: xmain.exe

main.o: main.c
	$(CC) -c $(OPTS) -I$(INCLUDECBLAS) main.c

xmain.exe: main.o
	$(LOADER) -o $@ main.o $(LIBCBLAS)

clean:
	rm *.o *.exe 
 
