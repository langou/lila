#include "lila.h"

int lila_dge_qr_ormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	double *tau;
	int info, jj;
	double *Aii, *Qij;
	int ml;
	double *V;
	double normV_square;
	int vl;
	
	Aii = A + i*lda + i;
	Qij = Q + j*ldq + i;
	ml = m - i;

	tau = (double *) malloc( k * sizeof(double));
//	for(jj = 0; jj < k; jj++) tau[jj] = T[jj+jj*ldt];

	for( jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++, vl--, V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', k, n, (0e+00), (0e+00), Qij, ldq );
	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'N', ml, n, k, Aii, lda, tau, Qij, ldq, work, lwork );

	free( tau );



/*

    work(1:k,1:n) = randn(k,n); % (we put randn to simulate an allocation as opposed to putting zeros)
    [ work ] = blas_gemm ( 'T', 'N', k, n, m-k-i+1, ( +1.e+00 ), A, i+k, i, lda, Q, i+k, j, ldq, ( 0.e+00 ), work, 1, 1, ldwork );
    [ work ] = blas_trmm ( 'L', 'U', 'N', 'N', k, n, ( +1.e+00 ), T, i, i, ldt, work, 1, 1, ldwork );
    Q( i:i+k-1,j:j+n-1 ) = work(1:k,1:n);
    [ Q ] = blas_trmm ( 'L', 'L', 'N', 'U', k, n, ( -1.e+00 ), A, i, i, lda, Q, i, j, ldq );
    [ Q ] = blas_gemm ( 'N', 'N', m-k-i+1, n, k, ( -1.e+00 ), A, i+k, i, lda, work, 1, 1, ldwork, ( +1.e+00 ), Q, i+k, j, ldq );


*/

	return 0;

}
