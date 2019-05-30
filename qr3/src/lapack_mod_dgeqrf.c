#include "qr2.h"

int lapack_mod_dgeqrf( int m, int n, int nb, double *A, int lda, double *tau, double *T, int ldt ){

	double *A11, *A12, *T11, *T12, *tau1;
	int ml, nl, ib, k, i;

	i = 0;
	k = n;

	ib = nb; if( k - i - nb < 0 ) ib = k - i;

	A11  = A;
	A12  = A+ib*lda;
	T11  = T;
	T12  = T+ib*ldt;
	tau1 = tau;

	ml = m;
	nl = n;

	while( k - i > nb ){

//		dgeqr2_( &ml, &ib, A11, &lda, tau1, work, &lwork );
//		lapack_ref_dgeqr2( ml, ib, A11, lda, tau1, work, lwork );
		lapack_mod_dgeqr2( ml, ib, A11, lda, tau1 );

		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, A11, lda, tau1, T11, ldt);

		LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'T', 'F', 'C', ml, nl-ib, ib, A11, lda, T11, ldt, A12, lda, T12, nl-ib);

		A11  += ib*(1+lda);
		A12  += ib*(1+lda);
		T11  += ib*ldt;
		T12  += ib*ldt;
		tau1 += ib;

		ml -= ib;		
		nl -= ib;		
		i  += ib;		

		ib = nb; if( k - i - nb < 0 ) ib = k - i;
	
	}	

//	dgeqr2_( &ml, &nl, A11, &lda, tau1, work, &lwork );
//	lapack_ref_dgeqr2( ml, nl, A11, lda, tau1, work, lwork );
	lapack_mod_dgeqr2( ml, nl, A11, lda, tau1 );


//	this larft is not in lapack_ref_orgqr
//	and it is not needed by lapack_mod_orgqr either
//	it is be useful when the first orgqr2 of orgqr uses ``T``
//	LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, nl, A11, lda, tau1, T11, ldt);

	return 0;

}
