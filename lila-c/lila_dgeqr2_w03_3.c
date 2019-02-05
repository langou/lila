#include "lila.h"

int lila_dgeqr2_w03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	//double *tau=NULL;
	double *work1;
	double *Aii, *Qii, *Tki;
	int j, k, ml;

	//tau = (double *) malloc( n * sizeof(double));
	work1 = (double *) malloc( n * n * sizeof(double));

	ml = m - i;
	k = i % mt;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tki = T + k + i*ldt;

  	info = dgeqr3( ml, n, Aii, lda, Tki, ldt );

// 	for(j = 0; j < n; j++) tau[j] = Tki[ j + j*ldt ];
//	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Qii, ldq);

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work1, n);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Tki, ldt, work1, n );

//	note that this is a triangle times triangle so we could have dtrtrmm() and dived # of FLOPS by 3x
	cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Aii, lda, work1, n );

//	note that the top n-by-n part is a lower times a upper and can be done dlauum
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (-1.0e+00), work1, n, Qii, ldq );

 	for(j = 0; j < n; j++) Qii[ j + ldq * j ] = 1.00e+00 + Qii[ j + ldq * j ];

	//free( tau );
	free( work1 );

	return 0;

}