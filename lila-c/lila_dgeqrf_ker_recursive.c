#include "lila.h"

int lila_dgeqrf_ker_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2, ml, vl, k, j; 
	double *tau=NULL, *Aii, *Qii, *V, normV_square;

	tau = (double *) malloc( n*n * sizeof(double));		

	Aii = A + i + i*lda;
	Qii = Q + i + i*ldq;

	ml = m - i;

	// This is a piece we'll need to make the good Q
	for( j = 0, vl=ml-1, V = Aii+1; j < k; j++, vl--, V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[j] = ( 2.0e+00 ) / normV_square ;
	}


	if ( n <= lila_param[ 3 ] ) {

  		info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );

	// This block below will work if  -   nx = mt = n     and everything else is commented out ...
	/*		
  		info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );
	*/

	} else {
	
		nb1 = n / 2;
		nb2 = n - nb1;

		info = lila_dgeqrf_ker_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		//lila_dormqrf
		info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', m, nb1, nb2, Aii, lda, tau, Qii, ldq, work, lwork );

		info = lila_dgeqrf_ker_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		
	
		//lila_dormqrbz     -     put the zeros above and then the magic
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', nb2, nb2, (0e+00), (0e+00), Qii, ldq );
		info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'N', ml, nb1, nb2, Aii, lda, tau, Qii, ldq, work, lwork );


	}
	
	free( tau );

	return 0;

}
