#include "qr3.h"

int our_dorgqr( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info ){

	double *Akk, *A0k, *Aik, *tauk;
	int kk, ml, nl, ib, i, j, flops, ldwork;
	
	ldwork = n;

	flops = 0;

	kk    = n;

	ml  = m-kk;
	nl  = n-kk;

	A0k  = A+kk*lda;
	Akk  = A+kk+kk*lda;
	tauk = tau+kk;

	ib = nb; if( kk - ib < 0 ) ib = n - nl;

	A0k  -= ib*lda;	
	Akk  -= ib*(1+lda);
	tauk -= ib;

	ml += ib;		
	nl += ib;		
	kk -= ib;

	for( i = 0; i < kk; i++){ for( j = 0; j < ib; j++ ){ A0k[i+j*lda] = (+0.0e00); } }

	flops += flops_org2r( ml, ib, ib );
	dorg2r_( &ml, &ib, &ib, Akk, &lda, tauk, work, &lwork );

	while( kk > 0 ){

		ib = nb; if( kk - ib < 0 ) ib = n - nl;

		A0k  -= ib*lda;
		Akk  -= ib*(1+lda);
		Aik   = Akk+ib*lda;
		tauk -= ib;

		ml += ib;		
		nl += ib;		
		kk -= ib;
	
		flops += flops_larft( ml, ib );
		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, Akk, lda, tauk, work, ldwork);

		flops += flops_larfb( ml, nl-ib, ib );
		info = LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', ml, nl-ib, ib, Akk, lda, work, ldwork, Akk+ib*lda, lda, work+ib, ldwork);

		flops += flops_org2r( ml, ib, ib );
		dorg2r_( &ml, &ib, &ib, Akk, &lda, tauk, work, &lwork );

		for( i = 0; i < kk; i++){ for( j = 0; j < ib; j++ ){ A0k[i+j*lda] = (+0.0e00); } }

	}	

	printf("flops in total = %d\n", flops);
	return flops;

}
