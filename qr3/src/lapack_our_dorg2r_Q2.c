#include "qr2.h"

int lapack_our_dorg2r_Q2( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *Axx, *tau1;
	int m1, nx, i, j;
	double tmp;

	for ( j = k; j < n; j++ ){
		for ( i = 0; i < m; i++ ){
			A[ i + j * lda ] = (+0.00e+00);
		}
		A[ j + j * lda ] = (+1.00e+00);
	}
		
	m1 = m-k;
	nx = n-k;

	A11 = A+(k)*(1+lda);
	tau1 = tau+k;
	Axx = A11;

	for ( j = k-1; j >=0; j-- ){

		A11  -= (1+lda);
		tau1 --;
		Axx --;
		m1 ++;		
	
		tmp = (*A11);
		(*A11) = (+1.00);
		qr2_aux_dlarf_wrapper( 'L', m1, nx, A11, 1, (*tau1), Axx, lda, work );
		(*A11) = tmp;

	}

	return 0;

}
