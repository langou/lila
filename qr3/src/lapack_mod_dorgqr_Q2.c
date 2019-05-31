#include "qr2.h"

int lapack_mod_dorgqr_Q2( int m, int n, int k, int nb, double *A, int lda, double *Q, int ldq, double *tau, double *work, int lwork ){

	double *A11, *Qxx, *Q0x, *tau1;
	int k0, m1, n1, nx, ib, ldwork;
	
	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	A11  = A+k0*(1+lda);
	tau1 = tau+k0;

	Qxx = Q+k0;
	nx  = n-k;

	LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

	{ int i, j; for( i = 0; i < k0; i++){ for( j = 0; j < nx; j++ ){ Q[i+j*ldq] = (+0.0e00); } } }

	lapack_mod_dorg2r_Q2( m1, n1, ib, A11, lda, work, ib, Qxx, ldq );

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A11  -= ib*(1+lda);
		Qxx  -= ib;
		tau1 -= ib;

		m1 += ib;		
		k0 -= ib;		
	
		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

		lapack_ref_dlarfb_lnfc( m1, nx, ib, A11, lda, work, ib, Qxx, lda, work+ib*ib );

//		lapack_mod_dorg2r_Q2( m1, nx, ib, A11, lda, work, ib, Qxx, ldq );
//		lapack_mod_dorg2r_Q2( m1, n1, ib, A11, lda, work, ib, Qxx, ldq );

	//	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, nx, m1-ib, (+1.0e+00), A11+ib, lda, Qxx+ib, ldq, (+0.0e+00), Qxx, ldq );
	//	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, ib, nx, (+1.0e+00), A11, lda, Qxx, ldq );
	//	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m1-ib, nx, ib, (-1.0e+00), A11+ib, lda, Qxx, ldq, (+1.0e+00), Qxx+ib, ldq );
	//	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, nx, (-1.0e+00), A11, lda, Qxx, ldq );


	}	

	return 0;

}
