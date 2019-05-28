#include "qr3.h"

int UinvLTinU_cheat( int n, double *L, int ldl, double *U, int ldu ){

	int i, j; 
	double *X;

	X = (double *) malloc( n * n * sizeof(double));

	for(i=0;i<n;i++){for(j=0;j<n;j++){ X[i+j*n] = L[j+i*ldl]; }}	
	LAPACKE_dlaset( LAPACK_COL_MAJOR, 'L', n, n, (0.0e+00), (1.0e+00), X, n);

	cblas_dtrsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, (+1.0e+00), U, ldu, X, n );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, X, n, U, ldu );

	free( X );

	return 0;

}
