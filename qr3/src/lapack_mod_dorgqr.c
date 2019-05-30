#include "qr2.h"

int lapack_mod_dorgqr( int m, int n, int k, int nb, double *A, int lda, double *tau, double *T, int ldt){

//	case 1: this routine can either use a sequence of T matrices provided by geqrf_mod (and passed through T)
//	case 2: it can creates its own T's one step at a time, in this second case T is in effect a nb-by-nb workspace, and ldt=nb, then do not move pointer in T!

	double *A11, *A01, *A12, *T1, *tau1;
	int k0, m1, n1, n2, ib;
//	int ldwork;
	
//	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	A01 = A+k0*lda;
	A11 = A+k0*(1+lda);
	T1  = T+k0*ldt;

	tau1 = tau+k0;

//	dorg2r_( &m1, &n1, &ib, A11, &lda, tau1, work, &lwork );
//	lapack_ref_dorg2r( m1, n1, ib, A11, lda, tau1, work, lwork );
	lapack_mod_dorg2r( m1, n1, ib, A11, lda, tau1 );

//	for( i = 0; i < k0; i++){ for( j = 0; j < n1; j++ ){ A01[i+j*lda] = (+0.0e00); } }

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A01  -= ib*lda;
		A11  -= ib*(1+lda);
		A12   = A11+ib*lda;
		T1   -= ib*ldt;
		tau1 -= ib;

		n2  = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		
	
//		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ldwork);
//		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, T1, ldt);

//		LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', m1, n2, ib, A11, lda, work, ldwork, A12, lda, work+ib, ldwork);
//		lapack_ref_dlarfb_lnfc( m1, n2, ib, A11, lda, work, ib, A12, lda, work+ib*ib );
		lapack_mod_dlarfb_lnfc_bz( m1, n2, ib, A11, lda, T1, ldt, A12, lda );

//		dorg2r_( &m1, &ib, &ib, A11, &lda, tau1, work, &lwork );
//		lapack_ref_dorg2r( m1, ib, ib, A11, lda, tau1, work, lwork );
		lapack_mod_dorg2r( m1, ib, ib, A11, lda, tau1 );

//		for( i = 0; i < k0; i++){ for( j = 0; j < ib; j++ ){ A01[i+j*lda] = (+0.0e00); } }

	}	

	return 0;

}
