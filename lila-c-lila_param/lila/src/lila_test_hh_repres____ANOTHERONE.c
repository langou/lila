#include "lila.h"

double lila_test_hh_repres( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *V, int ldv ){

	double *Aii, *Q, *Qii, *work;
	double orth_1, orth_2, repres_1, norm;
	int ml, info, lwork, jj, ii;

	check = 0.e+00;

	Q = (double *) malloc(m * (n+i) * sizeof(double));
	Aii = A+i+i*lda;
	Vii = V+i+i*ldv;
	Qii = Q+i+i*m;
	ml = m-i;

//	explicitly construct m-by-n matrix Q with orthonormal columns
	lila_param[1]=1; lila_param[2]=1; lila_param[4]=2;
	lila_dgeqrf( lila_param, m, n, i, mt, V, ldv, T, ldt, Q, m, work, lwork );

//	checking || I - Q^T * Q ||
	orth = lila_test_qq_orth_1( m, n, i, Q, m );
	check = LAPACKE_dlapy2( check, orth);

//	checking || A - Q * R || / || A ||
	norm_repres = lila_test_qr_repres_1( m, n, i, A, lda, Q, m, V, ldv );
	check = LAPACKE_dlapy2( check, norm_repres);

	free( Q );
	free( work );

	return check;
}
