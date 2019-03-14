#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, ii, lwork, ml, verbose, testing;
	double *A, *Q, *T, *As, *work=NULL, *Aii, *Qii, *Tii;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0);

    	m         = 87;
    	n         = 53;
	ii        = 6;
	lda       = -1;
	ldq       = -1;
	verbose   = 0;
	testing   = 0;

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
		if( strcmp( *(argv + i), "-ii") == 0) {
			ii  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-verbose") == 0) {
			verbose  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-testing") == 0) {
			testing  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( m < n+ii ){ printf("\n\n YOUR CHOICE OF n AND ii HAVE MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	ldt            = n+ii;

	if( verbose == 1 ){

		printf("m = %4d, ",         m);
		printf("ii = %4d, ",       ii);
		printf("n = %4d, ",         n);
		printf("lda = %4d, ",     lda);
		printf("ldq = %4d, ",     ldq);
	}

	A  = (double *) malloc(lda * (n+ii) * sizeof(double));
	As = (double *) malloc(lda * (n+ii) * sizeof(double));
	Q  = (double *) malloc(ldq * (n+ii) * sizeof(double));
	T  = (double *) malloc(ldt * (n+ii) * sizeof(double));

 	for(i = 0; i < lda * (n+ii); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * (n+ii); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	Tii   = T + ii + ii*ldt;
	Aii   = A + ii + ii*lda;
	Qii   = Q + ii + ii*ldq;
	ml    = m - ii;
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+ii, A, lda, As, lda );
	normA = LAPACKE_dlange_work   ( LAPACK_COL_MAJOR, 'F', ml, n, Aii, lda, work );

//	1
//	work = (double *) malloc( 1 * sizeof(double));
//	lwork = -1;
//	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );
//	lwork = (int) work[0]; 
//	free( work );
//	
//	tau = (double *) malloc( (n+ii) * sizeof(double));
//	work = (double *) malloc( lwork * sizeof(double));
//	
//	gettimeofday(&tp, NULL);
//	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
//	
//	info = dgeqr3( ml, n, Aii, lda, Tii, ldt );
//	for (i=0;i<n;i++) tau[i]=Tii[i+i*ldt];
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
//	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );
//	
//	gettimeofday(&tp, NULL);
//	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
//	
//	free( tau  );
//	free( work );

//	2
	lwork = n * n;
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	info = dgeqr3( ml, n, Aii, lda, Tii, ldt );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Qii, ldq);
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work, n);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Tii, ldt, work, n );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Aii, lda, work, n );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (-1.0e+00), work, n, Qii, ldq );
 	for(j = 0; j < n; j++) Qii[ j + ldq * j ] = 1.00e+00 + Qii[ j + ldq * j ];

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	free( work );

	if ( verbose == 0 ){

		printf("%6d %6d %16.8f %10.3f\n", m, n, elapsed_refL, perform_refL);
		//printf("%6d %6d %s %16.8f %10.3f\n", m, n, getenv("OPENBLAS_NUM_THREADS"), elapsed_refL, perform_refL);

	}
	if( testing == 1 ){

	double *QQ, *RR, *HH, norm_repres_1, norm_orth_1, *As_ii;
	double norm_orth_3, norm_repres_3, norm_diffQ_3, norm_repres_2_2, norm_orth_2, norm_repres_2_1;

	As_ii = As+ii+ii*lda;
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

	RR = (double *) malloc(m * n * sizeof(double));
	double *RRi0;
	RRi0 = RR+ii;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, As_ii, lda, RRi0, m );

	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
	if( m == n+ii ) lila_dormqrf_z03( m, n, n-1, ii, 0, n+ii, A, lda, T, ldt, RR, m, work, lwork );
	if ( m > n+ii ) lila_dormqrf_z03( m, n,   n, ii, 0, n+ii, A, lda, T, ldt, RR, m, work, lwork );
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
	if( m == n+ii ) lila_dormqrf_z00( ml, ml, n-1,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	if ( m > n+ii ) lila_dormqrf_z00( ml, ml,   n,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
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

	lwork = ml*n;
	work  = (double *) malloc(ml * n * sizeof(double));
 	for(i = 0; i < ml; i++) for(j = 0; j < n; j++) work[ i+j*ml ] = Qii[ i+j*ldq] - QQ[ i+j*ml ];
	norm_diffQ_3 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', ml, n, work, ml, NULL );
	free( work );
//	printf("9 \n");

	free( QQ );
	free( HH );


	printf("\n");
	printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);
	printf("\n");
	printf("| res1  = %5.1e    orth1 = %5.1e ", norm_repres_1, norm_orth_1);
	printf("\n");
	printf("| res2  = %5.1e    res2  = %5.1e  orth2 = %5.1e  ", norm_repres_2_1, norm_repres_2_2, norm_orth_2);
	printf("\n");
	printf("| res3  = %5.1e    orth3 = %5.1e  diff3 = %5.1e", norm_repres_3, norm_orth_3, norm_diffQ_3 );
	printf("\n");
	}

	free( Q );
	free( A );
	free( As );
	free( T );


	return 0;
}
