
include ../make.inc

CHECK = check_qq_orth.o      \
	check_qr_repres.o    \

all: $(CHECK)
	ar -cr libcheck.a $^
	ranlib libcheck.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *a

