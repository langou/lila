#include "qr2.h"

int lapack_mod_dorg2r( int m, int n, int k, double *A, int lda, double *tau ){

	double *A11, *A01, *tau1;
	int m1, n1, i, j;
	double *A12, *A22;
	double *A21;


	for ( j = k; j < n; j++ ){
		for ( i = 0; i < m; i++ ){
			A[ i + j * lda ] = (+0.00e+00);
		}
		A[ j + j * lda ] = (+1.00e+00);
	}
		
	m1 = m-k;
	n1 = n-k;

	A01 = A+(k)*lda;
	A11 = A+(k)*(1+lda);
	tau1 = tau+k;

	for ( j = k-1; j >=0; j-- ){

		A01  -= lda;
		A11  -= (1+lda);
		A12 = A11 + lda;
		A22 = A11 + lda + 1;
		A21 = A11 + 1;

		tau1 --;

		m1 ++;		
		n1 ++;
	
		if( 1 < n1 ) {

//			(*A11) = (+1.00);
//			qr2_aux_dlarf_wrapper( 'L', m1, n1-1, A11, 1, (*tau1), A11+lda, lda, work );

//			lapack_ref_dlarfb_lnfc( m1, n1-1, 1, A11, lda, tau1, 1, A11+lda, lda, work );

//			lapack_mod_dlarfb_lnfc_bz( m1, n1-1, 1, A11, lda, tau1, 1, A11+lda, lda );

////			A12 = A21^T * A22 
//			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, 1, n1-1, m1-1, (+1.0e+00), A21, lda, A22, lda, (+0.0e+00), A12, lda );
////			A12 = tau1 * A12
//			cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1, n1-1, (+1.0e+00), tau1, 1, A12, lda );
////			A22 = A22 - A21 * A12
//			cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m1-1, n1-1, 1, (-1.0e+00), A21, lda, A12, lda, (+1.0e+00), A22, lda );
////			A12 = - A12
//			cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1, n1-1, (-1.0e+00), A11, lda, A12, lda );
////			A12 = A21^T * A22 
//			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, 1, n1-1, m1-1, (+1.0e+00), A21, lda, A22, lda, (+0.0e+00), A12, lda );

//			A12^T = A22^T * A21
			cblas_dgemv( CblasColMajor, CblasTrans, m1-1, n1-1, (+1.0e+00), A22, lda, A21, 1, (+0.0e+00), A12, lda );

//			A12 = - tau1 * A12
			cblas_dscal( n1-1, -(*tau1), A12, lda );

//			A22 = A22 + A21 * A12
			cblas_dger( CblasColMajor, m1-1, n1-1, (+1.0e+00), A21, 1, A12, lda, A22, lda);

		}


		if( 1 < m1 ) cblas_dscal( m1-1, -(*tau1), A11+1, 1 );
     		(*A11) = (+1.0e+00) - (*tau1);

//		for( i = 0; i < j; i++){ A01[i] = (+0.0e00); }

	}	

	return 0;

}
