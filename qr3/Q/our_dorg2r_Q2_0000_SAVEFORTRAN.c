#include "qr3.h"

//extern void dlarf_( char *side, int *m, int *n, double *v, int *incv, double *tau, double *c, int *ldc, double *work);

int our_dorg2r_Q2( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork ){

	int i, j, l;

//
//	Initialise columns k+1:n to columns of the unit matrix
//

	for ( j = k; j < n; j++ ){
		for ( l = 0; l < m; l++ ){
			A[ l + j * lda ] = (+0.00e+00);
		}
		A[ j + j * lda ] = (+1.00e+00);
	}
		
	for ( i = k-1; 0 <= i; i-- ){
//
//		Apply H(i) to A(i:m,i:n) from the left
//
		if( i < n-1 ){

			A[ i + i*lda ] = (+1.0e+00);

			char charL;
			double *A11, *A12, *tau1;
			int m1, n1, ione;
			charL = 'L'; m1 = m-i; n1 = n-i-1; A11 =  &( A[ i + i*lda ] ); ione = 1; tau1= &(tau[i]);  A12 =  &( A[ i + (i+1)*lda ] );
		
//			dlarf_( &charL, &m1, &n1, A11, &ione, tau1, A12, &lda, work );

			our_dlarfb_lnfc( m1, n1, 1, A11, lda, tau1, 1, A12, lda, work );

		
		}

		if( i < m-1 )
			cblas_dscal( m-i-1, -tau[ i ], &( A[ i + 1 + i * lda ] ), 1 );
         	A[ i + i*lda ] = (+1.0e+00) - tau[ i ];

//
//		Set A(1:i-1,i) to zero
//
		for ( l = 0; l < i-1; l++ ){
			A[ l + i * lda ] = (+0.00e+00);
		}

	}

	return 0;

}
