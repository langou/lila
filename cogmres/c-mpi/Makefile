include Makefile.opts

%.o: %.c
	${CC} ${OPTS} ${CC_FLAGS} -c $*.c -o $@

SOURCES = test.c  \
		cgs_v0.c cgs_v1.c \
		mgs_v0.c mgs_v1.c \
		rowmgs_v0.c rowmgs_v1.c \
		choleskyqr_A_v0.c choleskyqr_A_v1.c choleskyqr_A_v2.c choleskyqr_A_v3.c \
		choleskyqr_B_v0.c choleskyqr_B_v1.c choleskyqr_B_v2.c choleskyqr_B_v3.c \
		reducehouseholder_A_v0.c \
		reducehouseholder_A_v1.c \
		reducehouseholder_A_v2.c \
		reducehouseholder_B_v0.c \
		reducehouseholder_B_v1.c \
		reducehouseholder_B_v2.c \
		reducehouseholder_B_v3.c \
		myqrfactorizationofTypeIIImatrices_A.c \
		myqrfactorizationofTypeIIImatrices_B.c \
		check_representativity.c check_orthogonality.c check_orthogonality_RFP.c check_Rfactor.c \
		LILA_mpiop_sum_upper.c \
		LILA_mpiop_qr_upper_v0.c LILA_mpiop_qr_upper_v1.c LILA_mpiop_qr_upper_v2.c \
		LILA_mpiop_dlapy2.c \
		LILA_qr_uppers.c \
		LILA_qr_dormqr.c \
		LILA_up.c  \
		LILA_down.c \
		dgeqr3.c
    

#		scalapackqr2_A.c scalapackqr2_B.c scalapackqrf_A.c scalapackqrf_B.c \

OBJECTS := $(SOURCES:%.c=%.o)

all: xtest

xtest: ${OBJECTS}
	${LOADER} ${OBJECTS} \
		${SCALAPACKLIB} \
		${LAPACKLIB} ${BLASLIB} ${LIBS} \
		-o $@
#
# clean
#
clean:
	rm -f x* *.o MATLAB_*.dat MATLAB_*.m
