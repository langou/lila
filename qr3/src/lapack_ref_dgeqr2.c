#include "qr2.h"

int lapack_ref_dgeqr2( int m, int n, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *A12, *tau1;
	int ml, nl;
	int i, j;
	double tmp;

	A11 = A;
	A12 = A+lda;
	tau1 = tau;

	ml = m;
	nl = n;

	while( nl > 1 ){

		LAPACKE_dlarfg_work( ml, A11, A11+1, 1, tau1 );

		tmp = (*A11);
		(*A11) = 1.0e+00;

		qr2_aux_dlarf_wrapper( 'L', ml, nl-1, A11, 1, (*tau1), A12, lda, work);

		(*A11) = tmp;

		A11 += (1+lda);
		A12 += (1+lda);
		tau1 ++;

		ml --;		
		nl --;		

	}	

	LAPACKE_dlarfg_work( ml, A11, A11+1, 1, tau1 );

	return 0;

}
