#include "qr3.h"

int our_dorg2r_BZ( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *A01, *tau1;
	int m1, n1, i, j;

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
		tau1 --;

		m1 ++;		
		n1 ++;
	
//		if( 1 < n1 ) {
//			(*A11) = (+1.00);
//			wrapper_dlarf( 'L', m1, n1-1, A11, 1, (*tau1), A11+lda, lda, work );
//		}

		
		//(*A11) = (+1.00);
		cblas_dgemv( CblasColMajor, CblasTrans, m1, n1-1, (*tau1), A11, lda, A11+lda, lda, (+1.0e00), work, lda);
		cblas_dger ( CblasColMajor, m1, n1-1, (+1.0e00), A11, lda, work, lda, A11, lda);






		if( 1 < m1 ) cblas_dscal( m1-1, -(*tau1), A11+1, 1 );
     		(*A11) = (+1.0e+00) - (*tau1);

		for( i = 0; i < j; i++){ A01[i] = (+0.0e00); }

	}	

	return 0;

}
