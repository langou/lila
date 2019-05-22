
include ../make.inc

SRC = check_qq_orth.o        \
	check_qr_repres.o    \

CC = gcc
OPTS = -Wall

all: $(SRC)
	ar -cr libcheck.a $^
	ranlib libcheck.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *a

