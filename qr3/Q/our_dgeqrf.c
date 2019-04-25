#include "qr3.h"

int our_dgeqrf( int m, int n, int nb, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *A12, *tau1;
	int ml, nl, ib, k, i, ldwork;
	
	ldwork = n;

	i = 0;
	k = n; // k = min(m,n)

	ib = nb; if( k - i - nb < 0 ) ib = k - i;

	A11  = A;
	A12  = A+ib*lda;
	tau1 = tau;

	ml = m;
	nl = n;

	while( k - i > nb ){

		dgeqr2_( &ml, &ib, A11, &lda, tau1, work, &lwork );

		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, A11, lda, tau1, work, ldwork);

		LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'T', 'F', 'C', ml, nl-ib, ib, A11, lda, work, ldwork, A12, lda, work+ib, ldwork);

		A11  += ib*(1+lda);
		A12  += ib*(1+lda);
		tau1 += ib;

		ml -= ib;		
		nl -= ib;		
		i  += ib;		

		ib = nb; if( k - i - nb < 0 ) ib = k - i;
	
	}	

	dgeqr2_( &ml, &nl, A11, &lda, tau1, work, &lwork );


	return 0;

}
