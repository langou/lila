#include <stdio.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

int lila_dge_qr_wq_vL0( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int j, info ; 
	double *tau=NULL;

	tau = (double *) malloc( n * sizeof(double));

	double *Aii, *Tii, *Qii;
	int ml;
	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + i*ldt + i;
	ml = m - i;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
  	for(j = 0; j < n; j++) Tii[j+j*ldt] = tau[j];

//	info = dgeqr3( ml, n, Aii, lda, Tii, ldt );
//	for(j = 0; j < n; j++) tau[j] = Tii[j+j*ldt];

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	free( tau );

	return 0;

}
