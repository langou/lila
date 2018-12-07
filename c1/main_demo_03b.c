#include "lila.h"

int main(int argc, char ** argv){

	int m, n, mt, nb, vb, i, j, l, lda, ldq, ldt, info, lwork;
	int *S;
	double *A, *As, *T, *Q, *work=NULL;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0); 

    	m = 30;
    	n = 20;
	nb = 5;
	mt = n;
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
		if( strcmp( *(argv + i), "-mt") == 0) {
			mt  = atoi( *(argv + i + 1) );
			i++;
		}

	}

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	printf("m = %4d; ",m);
	printf("n = %4d; ",n);
	printf("nb = %4d; ",nb);
	printf("mt = %4d; ",mt);
	printf("lda = %4d; ",lda);
	printf("ldq = %4d; ",ldq);
	printf("\n");

	A = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q = (double *) malloc(ldq * n * sizeof(double));
	S = (int *) malloc(n * sizeof(int));


 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

	ldt = mt;
	T = (double *) malloc(ldt * n * sizeof(double));

	lwork = 1920000;
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	// Cholesky QR on all of A
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq ); 
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, 0e+00, A, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, A, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1.0e+00, A, lda, Q, ldq );




	double *TTT, *AAA, *QQQ;
	int *SSS;
	TTT = (double *) malloc(n * n * sizeof(double));
	AAA = (double *) malloc(lda * n * sizeof(double));
	QQQ = (double *) malloc(ldq * n * sizeof(double));
	SSS = (int *) malloc(n * sizeof(int));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, AAA, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, QQQ, ldq ); 




	printf("\n");
	i = 0;
	j = 0;
	l = 0;
	if ( nb > n ) vb = n; else vb = nb;

		lila_dorghr_w03( m, vb, i, j, l, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

		lila_ormhr_w0b( m, vb, i, j, AAA, lda, TTT, n, QQQ, ldq, SSS );
		lila_dorgh2( m, vb, i, j, -1, AAA, lda, TTT, n, QQQ, ldq, work, lwork, SSS );
		info = lila_dlarft_connect_w02(m, vb, j, 0, -1, AAA, lda, TTT, n );

		j = vb;

	if ( j+nb > n ) vb = n-j; else vb = nb;
	while( vb!=0 ){

		lila_dorghr_w03( m, vb, i, j, l, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

		lila_ormhr_w0b( m, vb, i, j, AAA, lda, TTT, n, QQQ, ldq, SSS );
		lila_dorgh2( m, vb, i, j, -1, AAA, lda, TTT, n, QQQ, ldq, work, lwork, SSS );
		info = lila_dlarft_connect_w02(m, vb, j, 0, -1, AAA, lda, TTT, n );

		j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	}


//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, AAA, lda, A, lda ); 
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, QQQ, ldq, Q, ldq ); 

//printf("\nAAA - A =");
//for(i = 0; i < m; i++){ for(j = 0; j < n; j++){ AAA[ i + j*lda] = AAA[ i + j*lda] - A[ i + j*lda]; }  }
//printf("\n");
//for(i = 0; i < m; i++){ for(j = 0; j < n; j++){ printf(" %+5.1e ", AAA[ i + j*lda]); } printf("\n"); }
//printf("\n");

//printf("\nQQQ - Q = ");
//for(i = 0; i < m; i++){ for(j = 0; j < n; j++){ QQQ[ i + j*lda] = QQQ[ i + j*lda] - Q[ i + j*lda]; }  }
//printf("\n");
//for(i = 0; i < m; i++){ for(j = 0; j < n; j++){ printf(" %+5.1e ", QQQ[ i + j*lda]); } printf("\n"); }
//printf("\n");


printf("\n");
for(i = 0; i < mt; i++){ for(j = 0; j < n; j++){ printf(" %+5.1e ", T[ i + j*ldt]); } printf("\n"); }
printf("\n");

printf("\n");
for(i = 0; i < n; i++){ for(j = 0; j < n; j++){ printf(" %+5.1e ", TTT[ i + j*n]); } printf("\n"); }
printf("\n");

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;

////////////////////////////////////////////////////////////////////////////////////////////////////
	double norm_orth_1, norm_repres_1;
	double *QQ, *RR, *HH, norm_repres_2_1, norm_repres_2_2, norm_orth_2;
	double norm_orth_3, norm_repres_3, norm_diffQ_3;

//	Check the orthogonality of Q
	work = (double *) malloc(n * n * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
	norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );

//	Check A = QR
	work = (double *) malloc(m * n * sizeof(double));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, work, m );
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) work[i+j*m] -= As[i+j*lda];
	norm_repres_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	norm_repres_1 = norm_repres_1 / normA;
	free( work );

//	Check H^TA = R by looking at zeros below nxn-block 
	RR = (double *) malloc(m * n * sizeof(double));
	work = (double *) malloc(m * m * sizeof(double));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, As, lda, RR, m );
//	if( m >= n )  lila_dormqrf_z03( m, n, n, 0, 0, mt, A, lda, RR, m, T, ldt, work, lwork ); 
	if( m > n )  lila_dormqrf_z00( m, n, n, 0, 0, -1, A, lda, RR, m, T, ldt, work, lwork );
	if( m == n ) lila_dormqrf_z00( m, n, n-1, 0, 0, -1, A, lda, RR, m, T, ldt, work, lwork );
	free( work );

	if( m > n )  norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', m-1, n, RR+1, m, NULL );
	if( m == n ) norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', m-1, n-1, RR+1, m, NULL );

//	Check that R computed in A is the same generated from H^TA
	work = (double *) malloc(n * n * sizeof(double));
 	for(i = 0; i < n; i++) for(j = 0; j < n; j++) work[i+j*n] = A[i+j*lda] - RR[i+j*m];
	norm_repres_2_2 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', n, n, work, n, NULL );
	free( work );
	free( RR );

//	Applies V to identity to construct H
	HH = (double *) malloc(m * m * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', m, m, (0e+00), (1e+00), HH, m );
	work = (double *) malloc(m * m * sizeof(double));
	if( m >= n )  lila_dormqrf_z00( m, m, n, 0, 0, -1, A, lda, HH, m, T, ldt, work, lwork );
	//if( m == n ) lila_dormqrf_z00( m, m, n-1, 0, 0, -1, A, lda, HH, m, T, ldt, work, lwork );
	free( work );

//	Check H^TH = I
	work = (double *) malloc(m * m * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', m, m, (0e+00), (1e+00), work, m );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, m, m, 1.0e+00, HH, m, -1.0e+00, work, m );
	norm_orth_2 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, m, work, m, NULL );
	free( work );

//	Transpose H into QQ
	QQ = (double *) malloc(m * n * sizeof(double));
	for( i = 0; i < m; i++){ for( j = 0; j < n; j++){ QQ[i+j*m] = HH[j+i*m]; } } 

//	Check QQ^TQQ = I
	work = (double *) malloc(n * n * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, QQ, ldq, -1.0e+00, work, n );
	norm_orth_3 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );

//	Check A = QQR
	work = (double *) malloc(m * n * sizeof(double));
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, QQ, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, work, m );
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) work[i+j*m] -= As[i+j*lda];
	norm_repres_3 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	norm_repres_3 = norm_repres_3 / normA;
	free( work );

//	Checking that QQ = Q
	work = (double *) malloc(m * n * sizeof(double));
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) work[i+j*m] = Q[i+j*ldq] - QQ[i+j*m];
	norm_diffQ_3 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', m, n, work, m, NULL );
	free( work );
	free( QQ );
	free( HH );

////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);

	printf("\n");
	printf("| res1    = %5.1e    orth1   = %5.1e ", norm_repres_1, norm_orth_1);

	printf("\n");
	printf("| res2-1  = %5.1e    res2-2  = %5.1e  orth2 = %5.1e  ", norm_repres_2_1, norm_repres_2_2, norm_orth_2);

	printf("\n");
	printf("| res3    = %5.1e    orth3   = %5.1e  diff3 = %5.1e", norm_repres_3, norm_orth_3, norm_diffQ_3 );

	printf("\n");

	free( A );
	free( As );
	free( Q );
	free( T );
	free( S );

	free( AAA );
	free( QQQ );
	free( TTT );
	free( SSS );

	return 0;
}
