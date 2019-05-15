#include "qr3.h"

//extern void dlarf_( char *side, int *m, int *n, double *v, int *incv, double *tau, double *c, int *ldc, double *work);
extern int wrapper_dlarf( char side, int m, int n, double *V, int incv, double tau, double *C, int ldc, double *work);



int our_dorg2r_Q2( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *A01, *A12, *tau1;
	int k0, m1, n1, n2, ib, i, j, ldwork;
	int nb = 1;


	for ( j = k; j < n; j++ ){
		for ( i = 0; i < m; i++ ){
			A[ i + j * lda ] = (+0.00e+00);
		}
		A[ j + j * lda ] = (+1.00e+00);
	}
		
	
	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-1;
	printf(" k0=%d\n",k0);
	printf(" ib=%d\n",ib);

	m1 = m-k+1;
	n1 = n-k+1;

	A01 = A+(k-1)*lda;
	A11 = A+(k-1)*(1+lda);
	tau1 = tau+k0;

//	our_dlarfb_lnfc( m1, n1-1, 1, A11, lda, tau1, 1, A11+lda, lda, work );

	(*A11) = (+1.00);
	wrapper_dlarf( 'L', m1, n1-1, A11, 1, (*tau1), A11+lda, lda, work );

	cblas_dscal( m1-1, -(*tau1), A11+1, 1 );
    	(*A11) = (+1.0e+00) - (*tau1);

	for( i = 0; i < k0; i++){ A01[i] = (+0.0e00); }

	while( k0 > 0 ){
//	for ( i = k-2; 0 <= i; i-- ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A01  -= lda;
		A11  -= (1+lda);
		A12   = A11+lda;
		tau1 --;

		m1 ++;		
		n1 ++;
		k0 --;		
	
//		our_dlarfb_lnfc( m1, n1-1, 1, A11, lda, tau1, 1, A11+lda, lda, work );

		(*A11) = (+1.00);
		wrapper_dlarf( 'L', m1, n1-1, A11, 1, (*tau1), A11+lda, lda, work );


//		dorg2r_( &m1, &ib, &ib, A11, &lda, tau1, work, &lwork );

		if( 1 < m1 ) cblas_dscal( m1-1, -(*tau1), A11+1, 1 );
     		(*A11) = (+1.0e+00) - (*tau1);

		for( i = 0; i < k0; i++){ A01[i] = (+0.0e00); }

	}	

	return 0;

}

int wrapper_dlarf( char side, int m, int n, double *V, int incv, double tau, double *C, int ldc, double *work){

	dlarf_( &side, &m, &n, V, &incv, &tau, C, &ldc, work);


	return 0;
}
