
include make.inc

#all: xmain_lapack_ref.exe 
all: xmain_ls_qr3.exe 

xmain_lapack_ref.exe: main_lapack_ref.o
	$(LOADER) -o $@ main_lapack_ref.o $(LIB)

xmain_ls_qr3.exe: main_ls_qr3.o dgeqr3_R.o
	$(LOADER) -o $@ main_ls_qr3.o dgeqr3_R.o $(LIB)

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *exe 
