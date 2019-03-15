#include "lila.h"

int lila_main_test( int *lila_param, int m, int n, int ii, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *As, double normA, double elapsed_refL, double perform_refL ){


	double *QQ, *RR, *HH, norm_repres_2_1, norm_repres_2_2, norm_orth_2, *Qii, *Tii, *Aii, *As_ii;
	double norm_orth_3, norm_repres_3, norm_diffQ_3, norm_orth_1, norm_repres_1, *work;
	int ml, info, lwork, j, i;

	ml = m-ii;

	As_ii = As+ii+ii*lda;
	Aii   =  A+ii+ii*lda;
	Tii   =  T+ii+ii*ldt;

	if( ldq > 0 ){
	Qii   =  Q+ii+ii*ldq;
	lwork = n*n;
	work  = (double *) malloc(n * n * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, -1.0e+00, work, n );
	norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );
//	printf("1 \n");

	lwork = ml*n;
	work  = (double *) malloc(ml * n * sizeof(double));
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Qii, ldq, work, ml );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (1.0e+00), Aii, lda, work, ml );
 	for(i = 0; i < ml; i++) for(j = 0; j < n; j++) work[ i+j*ml ] -= As_ii[ i+j*lda ];
	norm_repres_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, n, work, ml, NULL );
	norm_repres_1 = norm_repres_1 / normA;
	free( work );
//	printf("2 \n");
	}

	RR = (double *) malloc(m * n * sizeof(double));
	double *RRi0;
	RRi0 = RR+ii;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, As_ii, lda, RRi0, m );

	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
	if ( ( lila_param[1] == 1 )&&( lila_param[2] == 2) ){
		lila_dormqrf_z03( m, n, n, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	} else if( lila_param[1] == 2 ) {
		lila_dormqrf_z03( m, n, n, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	} else {
		if( m == n+ii ) lila_dormqrf_z03( m, n, n-1, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
		if ( m > n+ii ) lila_dormqrf_z03( m, n,   n, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	}
	free( work );
//	printf("3 \n");

	norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', ml-1, n, RRi0+1, m, NULL );

	lwork = n*n;
	work  = (double *) malloc(n * n * sizeof(double));
 	for(i = 0; i < n; i++) for(j = 0; j < n; j++) work[ i+j*n ] = Aii[ i+j*lda ] - RRi0[ i+j*m ];
	norm_repres_2_2 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', n, n, work, n, NULL );
	free( work );
	free( RR );
//	printf("4 \n");

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

	free( QQ );
	free( HH );

	printf("\n");
	printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);
	if( ldq > 0 ){
	printf("\n");
	printf("| res1  = %5.1e    orth1 = %5.1e  diff3 = %5.1e", norm_repres_1, norm_orth_1,norm_diffQ_3);
	}
	printf("\n");
	printf("| res2  = %5.1e    res2  = %5.1e  orth2 = %5.1e  ", norm_repres_2_1, norm_repres_2_2, norm_orth_2);
	printf("\n");
	printf("| res3  = %5.1e    orth3 = %5.1e ", norm_repres_3, norm_orth_3 );
	printf("\n");


	return 0;
}

