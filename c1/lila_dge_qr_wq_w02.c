#include "lila.h"

int lila_dge_qr_wq_w02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int j, info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii;
	int ml;

	tau = (double *) malloc( n * sizeof(double));

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + i*ldt + i;

	ml = m - i;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tii, ldt);

	info = lila_dlarft_w03( m, n, i, mt, A, lda, TTT, llldddttt, tau);

//  	info = dgeqr3( ml, n, Aii, lda, Tii, ldt );
//	for(j = 0; j < n; j++) tau[j] = Tii[j+j*ldt];

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	free( tau );

	return 0;

}
