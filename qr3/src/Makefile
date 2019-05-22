
include ../make.inc

SRC = qr2_dgeqr3R.o \
	qr2_dgeqr3R_ISW.o

CC = gcc
OPTS = -Wall

all: $(SRC)
	ar -cr libqr2.a $^
	ranlib libqr2.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *a

