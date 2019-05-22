#include "lila.h"

double lila_test_vt_repres( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *V, int ldv ){

	double *Aii, *Q, *Qii, *Vii, *work;
	double orth, norm_repres, check;
	int ml, info, lwork;

	check = 0.e+00;

	lila_param[1]=1; lila_param[2]=1; lila_param[4]=2;
	lwork = lila_query_dgeqrf_w03_recursive( lila_param, m, n, i, mt, A, lda, T, ldt, V, ldv, work=NULL, -1 );
	
	work = (double *) malloc( lwork * sizeof(double));
	Q    = (double *) malloc( m * (n+i) * sizeof(double));
	Aii  = A+i+i*lda;
	Vii  = V+i+i*ldv;
	Qii  = Q+i+i*m;
	ml   = m-i;

//	explicitly construct m-by-n matrix Q with orthonormal columns
	lila_dgeqrf_recursive( lila_param, m, n, i, mt, V, ldv, T, ldt, Q, m, work, lwork );

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
