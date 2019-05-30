
include ../make.inc

SRC = qr2_dgeqr3_R.o                \
	qr2_dgeqr3_R_ISW.o          \
	qr2_dgeqr3_R_UT.o           \
	qr2_dgeqr3_R_UT_ISW.o       \
	qr2_dorgqr3.o               \
	qr2_dorgqr3_UT.o            \
	qr2_dorgqr3_VT2Q.o          \
	qr2_dorgqr3_VS2Q.o          \
	qr2_dlarft3.o               \
	qr2_dlarft3_ISW.o           \
	qr2_dlarft3_ISW_V2T.o       \
	qr2_dlarft3_ISW_V2T_UT.o    \
	qr2_dlarft3_UT.o            \
	qr2_dlarft3_UT_ISW.o        \
	qr2_aux_dV2N.o              \
	qr2_aux_dV2tau.o            \
	qr2_aux_dN2T.o              \
	qr2_aux_dVS2Q.o             \
	qr2_aux_dVT2Q.o             \
	qr2_aux_dmLUinA.o           \
	qr2_aux_dULTinU.o           \
	qr2_aux_dApUBTinA.o         \
	qr2_aux_dUinvLTinU.o        \
	lapack_ref_dgeqrf.o         \
	lapack_ref_dorg2r.o         \
	lapack_ref_dorgqr.o         \
	lapack_mod_dorgqr.o         \
	lapack_our_dorgqr_Q2.o      \
	lapack_our_dorg2r_Q2.o      \
	lapack_ref_dlarfb_lnfc.o    \
	lapack_mod_dlarfb_lnfc_bz.o \
	lapack_mod_dorg2r.o         \
	lapack_ref_dgeqr2.o         \
	lapack_mod_dgeqr2.o         \
	lapack_mod_dgeqrf.o         \
	qr2_aux_dlarf_wrapper.o

all: $(SRC)
	ar -cr libqr2.a $^
	ranlib libqr2.a

.c.o:
	$(CC) $(OPTS) -I$(INCLUDE) -c -o $@ $<

clean:
	rm -f *o *a

