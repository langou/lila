#include "lila.h"

int lila_dgeqrf_LAPACK_appendcols( int m, int k, int n, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	double *tau, *Qii, *Aii, *Qji, *Aji;
	int info, ml;
	
	tau = T;
	ml  = m - k;

	Qii = Q + k*ldq + k;
	Aii = A + k*lda + k;

	Qji = Q + k*ldq;
	Aji = A + k*lda;

	if ( k!= 0 ) info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', m, n, k, A, lda, tau, Aji, ldq, work, lwork );

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau+k, work, lwork ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau+k, work, lwork );

	if( k != 0 ) info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', k, n, (0e+00), (0e+00), Qji, ldq );
	if( k != 0 ) info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'N', m, n, k, A, lda, tau, Qji, ldq, work, lwork );

	return 0;

}
