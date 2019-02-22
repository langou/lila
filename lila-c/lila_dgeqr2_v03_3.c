#include "lila.h"

int lila_dgeqr2_v03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *work1, *Aii, *Tki;
	int j, k, ml, info; 

//	work1 = (double *) malloc( n * n * sizeof(double));
	work1 = work;	

	ml = m - i;
	k  =  i%mt;

	Aii = A + i*lda + i;
	Tki = T + i*ldt + k;

  	info = dgeqr3( ml, n, Aii, lda, Tki, ldt );

// not sure if below is needed

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work1, n);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Tki, ldt, work1, n );

//	note that this is a triangle times triangle so we could have dtrtrmm() and dived # of FLOPS by 3x
	cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Aii, lda, work1, n );

	return 0;

}
