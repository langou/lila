#include "lila.h"

double lila_test_r_represe_2( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *R, int ldr ){

	double norm_repres, *work, *RR, *RRi0, *Rii, *Aii, norm_repres_2;
	int ml, lwork, info, ii, jj;

	RR = (double *) malloc(m * n * sizeof(double));

	RRi0 = RR+i;
	Aii  = A+i+i*lda;
	Rii  = R+i+i*ldr;

	ml = m-i;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, RRi0, m );

	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
	if ( ( lila_param[1] == 1 )&&( lila_param[2] == 2) ){
		lila_dormqrf_z03( m, n, n, i, 0, mt, R, ldr, T, ldt, RR, m, work, lwork );
	} else if( lila_param[1] == 2 ) {
		lila_dormqrf_z03( m, n, n, i, 0, mt, R, ldr, T, ldt, RR, m, work, lwork );
	} else {
		if( m == n+i ) lila_dormqrf_z03( m, n, n-1, i, 0, mt, R, ldr, T, ldt, RR, m, work, lwork );
		if ( m > n+i ) lila_dormqrf_z03( m, n,   n, i, 0, mt, R, ldr, T, ldt, RR, m, work, lwork );
	}
	free( work );
	norm_repres = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', ml-1, n, RRi0+1, m, NULL );

	lwork = n*n;
	work  = (double *) malloc(n * n * sizeof(double));
 	for(ii = 0; ii < n; ii++) for(jj = 0; jj < n; jj++) work[ ii+jj*n ] = Rii[ ii+jj*ldr ] - RRi0[ ii+jj*m ];
	norm_repres_2 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', n, n, work, n, NULL );
	free( work );
	free( RR );


	return norm_repres, norm_repres_2;
}
