#include "qr2.h"

int lapack_mod_dgeqr2( int m, int n, double *A, int lda, double *tau ){

//	two modifications:
//	(1) we not do V11 =1.0e+00, this saves on scale, and enable better thread safety (well maybe not here)
//	(2) we set workspace to be inside tau (in the tau space that we will use later on)
//
//	tau needs to be of size n, in the case m>n
//	tau needs to be of size n-1, in the case m=n, maybe we want tau[n-1] = 0 in this case, not sure

	double *A11, *A21, *A12, *A22, *tau1, *W;
	int ml, nl;
	int j;

	A11 = A;
	A21 = A11 + 1;
	A12 = A+lda;
	A22 = A12 + 1;
	tau1 = tau;
	W = tau1 + 1;

	ml = m;
	nl = n;

	while( nl > 1 ){

		LAPACKE_dlarfg_work( ml, A11, A11+1, 1, tau1 );

		cblas_dcopy( nl-1, A12, lda, W, 1);
		cblas_dgemv( CblasColMajor, CblasTrans, ml-1, nl-1, (+1.0e+00), A22, lda, A21, 1, (+1.0e+00), W, 1);
		cblas_dscal( nl-1, (*tau1), W, 1);
		cblas_dger( CblasColMajor, ml-1, nl-1, (-1.0e+00), A21, 1, W, 1, A22, lda);
		for(j = 0; j < nl-1; j++) A12[ j*lda ] -= W[ j ];

		A11 += (1+lda);
		A21 = A11 + 1;
		A12 += (1+lda);
		A22 = A12 + 1;
		tau1 ++;
		W = tau1 + 1;

		ml --;		
		nl --;		

	}	

	LAPACKE_dlarfg_work( ml, A11, A11+1, 1, tau1 );

	return 0;

}
