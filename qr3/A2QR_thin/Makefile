
include ../make.inc

#all: xmain_A2QR_lapack.exe 
#all: xmain_A2QR_qr3.exe  
#all: xmain_A2QR_qr3_ISW.exe  
#all: xmain_A2QR_qr3_VT2Q.exe  
#all: xmain_A2QR_qr3_VT2Q_ISW.exe  

all: xmain_A2QR_lapack.exe xmain_A2QR_qr3.exe xmain_A2QR_qr3_ISW.exe xmain_A2QR_qr3_VT2Q.exe xmain_A2QR_qr3_VT2Q_ISW.exe  

xmain_A2QR_lapack.exe: main_A2QR_lapack.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_lapack.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3.exe: main_A2QR_qr3.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_ISW.exe: main_A2QR_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_VT2Q.exe: main_A2QR_qr3_VT2Q.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_VT2Q.o ../src/libqr2.a ../check/libcheck.a $(LIB)

xmain_A2QR_qr3_VT2Q_ISW.exe: main_A2QR_qr3_VT2Q_ISW.o ../src/libqr2.a ../check/libcheck.a
	$(LOADER) -o $@ main_A2QR_qr3_VT2Q_ISW.o ../src/libqr2.a ../check/libcheck.a $(LIB)

../src/libqr2.a:  
	make -C ../src/

../check/libcheck.a: 
	make -C ../check/

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe ../src/*a ../src/*o ../check/*a ../check/*o
