#include "lila.h"

int main(int argc, char ** argv){

	int m, n, nb, vb, mt, i, j, lda, ldq, ldt, info, lwork, lwork1;
	double *A, *As, *T, *Q, *work=NULL;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0);

    	m = 30;
    	n = 20;
	nb = 1;
	lda = -1;
	ldq = -1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	printf("m = %4d; ",m);
	printf("n = %4d; ",n);
	printf("nb = %4d; ",nb);
	printf("lda = %4d; ",lda);
	printf("ldq = %4d; ",ldq);
	printf("  ");

	A = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q = (double *) malloc(ldq * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

	ldt = n;
	T = (double *) malloc(ldt * n * sizeof(double));
	mt = -1;

	printf("\n");
	lwork = 0;
	j = 0;
	if ( nb > n ) vb = n; else vb = nb;
	lwork1 =  lila_wsq_dgeqr2_w02a( m, vb, 0, mt, A, lda, T, ldt, Q, ldq, work ); if ( lwork1 > lwork ) lwork = lwork1; printf("GEQR2   lwork1 = %4d;\n",lwork1);
	j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	while( vb!=0 ){
	lwork1 =  lila_wsq_dormqrf_z02( m, vb, j, 0, j, mt, A, lda, T, ldt, A, lda, work ); if ( lwork1 > lwork ) lwork = lwork1; printf("ORMQRF  lwork1 = %4d;\n",lwork1);
	lwork1 =  lila_wsq_dgeqr2_w02a( m, vb, j,       mt, A, lda, T, ldt, Q, ldq, work ); if ( lwork1 > lwork ) lwork = lwork1; printf("GEQR2   lwork1 = %4d;\n",lwork1);
	lwork1 = lila_wsq_dormqrbz_w02( m, vb, j, 0, j, mt, A, lda, T, ldt, Q, ldq, work ); if ( lwork1 > lwork ) lwork = lwork1; printf("ORMQRBZ lwork1 = %4d;\n",lwork1);
	j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	}

	printf("lwork = %4d; ",lwork);
	
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	j = 0;
	if ( nb > n ) vb = n; else vb = nb;
	info =        lila_dgeqr2_w02a( m, vb, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	while( vb!=0 ){
	info =        lila_dormqrf_z02( m, vb, j, 0, j, mt, A, lda, T, ldt, A, lda, work, lwork );
	info =        lila_dgeqr2_w02a( m, vb, j,       mt, A, lda, T, ldt, Q, ldq, work, lwork );
	info = lila_dlarft_connect_w02( m, vb, j, 0,    mt, A, lda, T, ldt );
	info =       lila_dormqrbz_w02( m, vb, j, 0, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );	
	j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	}

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;

////////////////////////////////////////////////////////////////////////////////////////////////////
	double norm_orth_1, norm_repres_1;
	double *QQ, *RR, *HH, norm_repres_2_1, norm_repres_2_2, norm_orth_2;
	double norm_orth_3, norm_repres_3, norm_diffQ_3;

	lwork = n*n;
	work = (double *) malloc(n * n * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
	norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );
	printf("1 \n");

	lwork = m*n;
	work = (double *) malloc(m * n * sizeof(double));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, work, m );
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) work[i+j*m] -= As[i+j*lda];
	norm_repres_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	norm_repres_1 = norm_repres_1 / normA;
	free( work );
	printf("2 \n");

	RR = (double *) malloc(m * n * sizeof(double));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, As, lda, RR, m );

	lwork = m*m;
	work = (double *) malloc(m * m * sizeof(double));
	if( m > n )  lila_dormqrf_z02( m, n, n, 0, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	if( m == n ) lila_dormqrf_z02( m, n, n-1, 0, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	free( work );
	printf("3 \n");

	if( m > n )  norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', m-1, n, RR+1, m, NULL );
	if( m == n ) norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', m-1, n-1, RR+1, m, NULL );

	lwork = n*n;
	work = (double *) malloc(n * n * sizeof(double));
 	for(i = 0; i < n; i++) for(j = 0; j < n; j++) work[i+j*n] = A[i+j*lda] - RR[i+j*m];
	norm_repres_2_2 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', n, n, work, n, NULL );
	free( work );

	free( RR );
	printf("4 \n");

	HH = (double *) malloc(m * m * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', m, m, (0e+00), (1e+00), HH, m );
	lwork = m*m;
	work = (double *) malloc(m * m * sizeof(double));
	if( m > n )  lila_dormqrf_z00( m, m, n, 0, 0, -1, A, lda, HH, m, T, ldt, work, lwork );
	if( m == n ) lila_dormqrf_z00( m, m, n-1, 0, 0, -1, A, lda, HH, m, T, ldt, work, lwork );
	free( work );
	printf("5 \n");

	lwork = m*m;
	work = (double *) malloc(m * m * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', m, m, (0e+00), (1e+00), work, m );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, m, m, 1.0e+00, HH, m, -1.0e+00, work, m );
	norm_orth_2 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, m, work, m, NULL );
	free( work );
	printf("6 \n");

	QQ = (double *) malloc(m * n * sizeof(double));

	for( i = 0; i < m; i++){ for( j = 0; j < n; j++){ QQ[i+j*m] = HH[j+i*m]; } } 

	lwork = n*n;
	work = (double *) malloc(n * n * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, QQ, ldq, -1.0e+00, work, n );
	norm_orth_3 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );
	printf("7 \n");

	lwork = m*n;
	work = (double *) malloc(m * n * sizeof(double));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, QQ, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, work, m );
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) work[i+j*m] -= As[i+j*lda];
	norm_repres_3 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	norm_repres_3 = norm_repres_3 / normA;
	free( work );
	printf("8 \n");

	lwork = m*n;
	work = (double *) malloc(m * n * sizeof(double));
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) work[i+j*m] = Q[i+j*lda] - QQ[i+j*m];
	norm_diffQ_3 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', m, n, work, m, NULL );
	free( work );
	printf("9 \n");

	free( QQ );

	free( HH );

////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);

	printf("\n");
	printf("| res1  = %5.1e    orth1 = %5.1e ", norm_repres_1, norm_orth_1);

	printf("\n");
	printf("| res2  = %5.1e    res2  = %5.1e  orth2 = %5.1e  ", norm_repres_2_1, norm_repres_2_2, norm_orth_2);

	printf("\n");
	printf("| res3  = %5.1e    orth3 = %5.1e  diff3 = %5.1e", norm_repres_3, norm_orth_3, norm_diffQ_3 );

	printf("\n");


	free( Q );
	free( A );
	free( As );
	free( T );

	return 0;
}
