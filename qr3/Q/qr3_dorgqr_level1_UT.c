#include "qr3.h"

int qr3_dorgqr_level1_UT( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *A01, *A12, *tau1;
	int k0, m1, n1, n2, ib;
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	A01 = A+k0*lda;
	A11 = A+k0*(1+lda);
	tau1 = tau+k0;

	dorg2r_( &m1, &n1, &ib, A11, &lda, tau1, work, &lwork );
//	dVS2Q( m1, ib, A11, lda  );

//	{ int i,j; for( i = 0; i < k0; i++){ for( j = 0; j < n1; j++ ){ A01[i+j*lda] = (+0.0e00); } } }

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A01  -= ib*lda;
		A11  -= ib*(1+lda);
		A12   = A11+ib*lda;
		tau1 -= ib;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		
	
////////////////

//		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, A11, lda);

		{ int i,j; for(i=0;i<ib;i++){ for(j=0;j<i;j++){ A11[j+i*lda] = A11[i+j*lda];}} }
		dV2N( ib, A11, lda );
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ib, m1-ib, (+1.0e+00), A11+ib, lda, (+1.0e+00), A11, lda );
		{ int i; for(i=0;i<ib;i++){ A11[i+i*lda] = 1/tau1[i]; } }



////////////////


////////////////

//		our_dlarfb_lnfc( m1, n2, ib, A11, lda, A11, lda, A12, lda, work );

		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, n2, m1-ib, (+1.0e+00), A11+ib, lda, A12+ib, lda, (+0.0e+00), A12, lda );
		cblas_dtrsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, ib, n2, (+1.0e+00), A11, lda, A12, lda );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m1-ib, n2, ib, (-1.0e+00), A11+ib, lda, A12, lda, (+1.0e+00), A12+ib, lda );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, n2, (-1.0e+00), A11, lda, A12, lda );

////////////////

////////////////

//		dorg2r_( &m1, &ib, &ib, A11, &lda, tau1, work, &lwork );

		dVS2Q( m1, ib, A11, lda  );

////////////////

//		{ int i,j; for( i = 0; i < k0; i++){ for( j = 0; j < ib; j++ ){ A01[i+j*lda] = (+0.0e00); } } }

	}	

	return 0;

}
