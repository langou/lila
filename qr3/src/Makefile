
include ../make.inc

SRC = qr2_dgeqr3R.o            \
	qr2_dgeqr3R_ISW.o      \
	qr2_dgeqr3R_UT.o       \
	qr2_dgeqr3R_UT_ISW.o   \
	qr2_dorgqr.o           \
	qr2_dorgqr_UT.o        \
	qr2_dorgqr_VT2Q.o      \
	qr2_dorgqr_VS2Q.o      \
	qr2_dV2tau.o           \
	qr2_larft.o            \
	qr2_larft_ISW.o        \
	qr2_larft_ISW_V2T.o    \
	qr2_larft_ISW_V2T_UT.o \
	qr2_larft_UT.o         \
	qr2_larft_UT_ISW.o     \
	qr2_dV2N.o             \
	qr2_dN2T.o             \
	qr2_dVS2Q.o            \
	qr2_dVT2Q.o            \
	qr2_mLUinA.o           \
	qr2_ULTinU.o           \
	qr2_ApUBTinA.o         \
	qr2_UinvLTinU.o        \

CC = gcc
OPTS = -Wall

all: $(SRC)
	ar -cr libqr2.a $^
	ranlib libqr2.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *a

