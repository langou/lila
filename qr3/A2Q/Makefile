
include ../make.inc

all: xmain_A2Q.exe xmain_A2Q_ISW.exe 

xmain_A2Q.exe: main_A2Q.o dgeqr3_Q.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_A2Q.o dgeqr3_Q.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

xmain_A2Q_ISW.exe: main_A2Q_ISW.o dgeqr3_Q.o dgeqr3_Q_ISW.o
	make -C ../src/
	make -C ../flops/
	make -C ../check/
	$(LOADER) -o $@ main_A2Q_ISW.o dgeqr3_Q.o  dgeqr3_Q_ISW.o ../src/libqr2.a ../check/libcheck.a ../flops/libflops.a $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o
