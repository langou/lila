#include "qr3.h"

int dorgqr( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info ){

	double *Akk, *T, *tauk;
	int kk, ldt, ml, nl, ib, i, j;
	
	// cheating
	ldt = nb;
	T   = (double *) malloc( ldt * nb * sizeof(double)); // note, I never move in T, I always overwrite my nb block

	kk  = n;

	ml  = m-kk;
	nl  = n-kk; // 0

	Akk  = A+kk+kk*lda;
	tauk = tau+kk;

	ib = nb; if( kk - ib < 0 ) ib = n - nl;
	
	Akk  -= ib*(1+lda);
	tauk -= ib;

	ml += ib;		
	nl += ib;		
	kk -= ib;
//
//	use unblocked code for the last or only block
//
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, ib, ib, Akk, lda, tauk, work, lwork );
//
//	Use blocked code
//
	while( kk > 0 ){

		ib = nb; if( kk - ib < 0 ) ib = n - nl;

		Akk  -= ib*(1+lda);
		tauk -= ib;

		ml += ib;		
		nl += ib;		
		kk -= ib;
	
//	
//		LAPACK method for constructing T
//
	//	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, Akk, lda, tauk, T, ldt);
//
//		Our 'faster' method for constructing T
//
		for(i=0;i<ib;i++){T[i+i*ldt] = 1.0e+00; for(j=0;j<i;j++){ T[j+i*ldt] = Akk[i+j*lda];}}
		xV2N( ib, T, ldt );
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ib, ml-ib, (+1.0e+00), Akk+ib, lda, (+1.0e+00), T, ldt );
		xN2T( ib, tauk, T, ldt );


		info = LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', ml, nl-ib, ib, Akk, lda, T, ldt, Akk+ib*lda, lda, work, lwork);

		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, ib, ib, Akk, lda, tauk, work, lwork );

	}	

	free( T );

	
//		This doesn't need to be in the while, but also not sure if it needs to be within this routine

		//if (m==n) T[n-1+ldt*(n-1)]=0.0e+00;


	return 0;

}
