#include "qr3.h"

int qr3_dA2QRTV_fake( int m, int n, double *A, int lda, double *Q, int ldq, double *T, int ldt ){

	double *tau;
	int i;

	tau = (double *) malloc( n * sizeof(double));

	dgeqr3( m, n, A, lda, T, ldt );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m, n, A, lda, Q, ldq );
	for(i=0;i<n;i++){tau[i] = T[i+i*ldt];}
	qr3_dorgqr( m, n, Q, ldq, T, ldt, tau );
	
	free( tau );

	return 0;

}
