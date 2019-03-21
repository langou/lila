#include "lila.h"

double lila_test_r_repres_2( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *R, int ldr ){

	double *HH;
	int ml, info;

	ml = m-i;

	HH = (double *) malloc(ml * ml * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', ml, ml, (0e+00), (1e+00), HH, ml );
	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
	if( ( lila_param[1] == 1 )&&( lila_param[2] == 2) ){
		lila_dormqrf_z00( ml, ml, n,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	} else if ( lila_param[1] == 2 ) {
		lila_dormqrf_z00( ml, ml, n,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	} else {
		if( m == n+ii ) lila_dormqrf_z00( ml, ml, n-1,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
		if ( m > n+ii ) lila_dormqrf_z00( ml, ml,   n,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	}
	free( work );
//	printf("5 \n");

	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', ml, ml, (0e+00), (1e+00), work, ml );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ml, ml, 1.0e+00, HH, ml, -1.0e+00, work, ml );
	norm_orth_2 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, ml, work, ml, NULL );
	free( work );
//	printf("6 \n");

	QQ = (double *) malloc(ml * n * sizeof(double));

	for( i = 0; i < ml; i++){ for( j = 0; j < n; j++){ QQ[i+j*ml] = HH[j+i*ml]; } } 

	lwork = n*n;
	work  = (double *) malloc(n * n * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, QQ, ml, -1.0e+00, work, n );
	norm_orth_3 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );
//	printf("7 \n");

	lwork = ml*n;
	work  = (double *) malloc(ml * n * sizeof(double));
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, QQ, ml, work, ml );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (1.0e+00), Aii, lda, work, ml );
 	for(i = 0; i < ml; i++) for(j = 0; j < n; j++) work[ i+j*ml ] -= As_ii[ i+j*lda];
	norm_repres_3 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, n, work, ml, NULL );
	norm_repres_3 = norm_repres_3 / normA;
	free( work );
//	printf("8 \n");

	if( ldq > 0 ){
	lwork = ml*n;
	work  = (double *) malloc(ml * n * sizeof(double));
 	for(i = 0; i < ml; i++) for(j = 0; j < n; j++) work[ i+j*ml ] = Qii[ i+j*ldq] - QQ[ i+j*ml ];
	norm_diffQ_3 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', ml, n, work, ml, NULL );
	free( work );
//	printf("9 \n");
	}





	return ;
}
