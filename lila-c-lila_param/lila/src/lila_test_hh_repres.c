#include "lila.h"

double lila_test_hh_repres( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *V, int ldv ){

//	this check creates and works on the full m-by-m matrix H
//	it is not recommended if m >> n

	double *H, *Hii, *Vii, *Aii, *work;
	double norm_orth, norm_repres, check;
	int ml, info, lwork, jj, ii;

	if( mt > n ) mt = n;
	check = 0.e+00;

	H = (double *) malloc(m * m * sizeof(double));

	Aii = A+i+i*lda;
	Vii = V+i+i*ldv;
	Hii = H+i+i*m;
	ml = m-i;

//	explicitly construct m-by-m unitary matrix H
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', ml, ml, (0e+00), (1e+00), Hii, m );
	lwork = ml * n;
//	work = mt * n; // seg fault for ./xmain_recursive.exe  -vrtq 0 -testing 1 -verbose 1 -m 5000 -n 2790 -nx 100 -mt 500
	work  = (double *) malloc( lwork * sizeof(double));
	if( ( lila_param[1] == 1 )&&( lila_param[2] == 2) ){
		lila_dormqrf_z03( m, m, n, i, 0, mt, V, ldv, T, ldt, H, m, work, lwork );
	} else if ( lila_param[1] == 2 ) {
		lila_dormqrf_z03( m, m, n, i, 0, mt, V, ldv, T, ldt, H, m, work, lwork );
	} else {
		if( m == n+i ) lila_dormqrf_z03( m, m, n-1,  i, 0, mt, V, ldv, T, ldt, H, m, work, lwork );
		if ( m > n+i ) lila_dormqrf_z03( m, m,   n,  i, 0, mt, V, ldv, T, ldt, H, m, work, lwork );
	}
	free(work);

//	checking || I - H^T * H ||
	work  = (double *) malloc( ml * ml * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', ml, ml, (0e+00), (1e+00), work, ml );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ml, ml, 1.0e+00, Hii, m, -1.0e+00, work, ml );
	norm_orth = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, ml, work, ml, NULL );
	check = LAPACKE_dlapy2( check, norm_orth);
	free(work);

//	checking || H^T A - [ R; 0 ] || / || A ||

	work  = (double *) malloc(ml * n * sizeof(double));
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Vii, ldv, work, ml );
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'L', ml-1, n, (0e+00), (0e+00), work+1, ml );

	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml, n, ml, -1.0e+00, Hii, m, Aii, lda, +1.0e+00, work, ml);

	norm_repres = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, n, work, ml, NULL );
	check = LAPACKE_dlapy2( check, norm_repres);

	free( work );
	free( H );

	return check;
}
