#include "null.h"

//int lapack_mod_dorgq2r_Q2( int m, int n, int k, int nb, double *A, int lda, double *Q, int ldq, double *tau, double *work, int lwork ){

int lapack_our_dorgq2r( int m, int n, int k, int nb, double *A, int lda, double *Q, int ldq, double *tau, double *work, int lwork ){

	double *A11, *Q1, *tau1;
	int k0, m1, n1, n2, ib, ldwork;
	
	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;
	n2 = n-k;

	A11  = A+k0*(1+lda);
	Q1   = Q+k0;
	tau1 = tau+k0;

	LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

//	lapack_mod_dorgq2r_Q2( m1, n1, ib, A11, lda, work, ib, Q1, ldq );
	qr3_aux_dorgq2r( m1, n1, ib, A11, lda, work, ib, Q1, ldq );

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A11  -= ib*(1+lda);
		Q1   -= ib;
		tau1 -= ib;

		m1 += ib;		
		k0 -= ib;		
	
		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, n2, m1-ib, (+1.0e+00), A11+ib, lda, Q1+ib, ldq, (+0.0e+00), Q1, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, ib, n2, (+1.0e+00), A11, lda, Q1, ldq );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m1-ib, n2, ib, (-1.0e+00), A11+ib, lda, Q1, ldq, (+1.0e+00), Q1+ib, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, n2, (-1.0e+00), A11, lda, Q1, ldq );

	}	

	return 0;

}
