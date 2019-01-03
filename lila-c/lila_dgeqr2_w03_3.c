#include "lila.h"

int lila_dgeqr2_w03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Qii, *Tki;
	int j, k, ml;

	tau = (double *) malloc( n * sizeof(double));

	ml = m - i;
	k = i % mt;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tki = T + k + i*ldt;

  	info = dgeqr3( ml, n, Aii, lda, Tki, ldt );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
 	for(j = 0; j < n; j++) tau[j] = Tki[ j + j*ldt ];
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	free( tau );

	return 0;

}
