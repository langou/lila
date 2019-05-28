
include ../make.inc

SRC = qr2_dgeqr3_R.o                \
	qr2_dgeqr3_R_ISW.o          \
	qr2_dgeqr3_R_UT.o           \
	qr2_dgeqr3_R_UT_ISW.o       \
	qr2_dorgqr3.o               \
	qr2_dorgqr3_UT.o            \
	qr2_dorgqr3_VT2Q.o          \
	qr2_dorgqr3_VS2Q.o          \
	qr2_dV2tau.o                \
	qr2_larft3.o                \
	qr2_larft3_ISW.o            \
	qr2_larft3_ISW_V2T.o        \
	qr2_larft3_ISW_V2T_UT.o     \
	qr2_larft3_UT.o             \
	qr2_larft3_UT_ISW.o         \
	qr2_dV2N.o                  \
	qr2_dN2T.o                  \
	qr2_dVS2Q.o                 \
	qr2_dVT2Q.o                 \
	qr2_mLUinA.o                \
	qr2_ULTinU.o                \
	qr2_ApUBTinA.o              \
	qr2_UinvLTinU.o             \
	level1_lapack_dorgqr_Q2.o   \
	our_dorg2r_Q2.o             \
	our_dlarfb_lnfc.o           \
	wrapper_dlarf.o

all: $(SRC)
	ar -cr libqr2.a $^
	ranlib libqr2.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *a

