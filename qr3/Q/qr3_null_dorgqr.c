#include "qr3.h"

int qr3_null_dorgqr( int m, int n, int k, int nb, double *A, int lda, double *Q, int ldq, double *tau, double *work, int lwork ){

	double *A11, *Q2, *tau1;
	int k0, m1, n1, n2, ib, i, j, ldwork;
	
	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;
	n2 = n-k;

	A11 = A+k0*(1+lda);
	Q2 = Q+k0;
	tau1 = tau+k0;

	LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);
	dorgqr_after( m1, n1, ib, A11, lda, work, ib, Q2, ldq );

//	for( i = 0; i < k0; i++){ for( j = 0; j < n2; j++ ){ Q[i+j*ldq] = (+0.0e00); } }

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A11  -= ib*(1+lda);
		Q2   -= ib;
		tau1 -= ib;

		m1 += ib;		
		k0 -= ib;		
	
		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

//		our_dlarfb_lnfc( m1, n2, ib, A11, lda, work, ib, Q2, ldq, work+ib*ib );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, n2, m1-ib, (+1.0e+00), A11+ib, lda, Q2+ib, ldq, (+0.0e+00), Q2, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, ib, n2, (+1.0e+00), A11, lda, Q2, ldq );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m1-ib, n2, ib, (-1.0e+00), A11+ib, lda, Q2, ldq, (+1.0e+00), Q2+ib, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, n2, (-1.0e+00), A11, lda, Q2, ldq );

	}	

	return 0;

}
