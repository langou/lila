#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, w, *S;
	double *A, *Q, *As, *T, *work=NULL;
	double normA, normR;
	double elapsed_refL, perform_refL;
	struct timeval tp;
	int n_lvl;
	int *nb_lvl;
	char mode;
	double norm_orth;

	srand(0);

    	m = 31;
    	n = 15;
	lda = -1;
	ldq = -1;

	mt = 9;
	n_lvl = 1;
	nb_lvl = (int *) malloc(n_lvl * sizeof(int));
	nb_lvl[0] = 10;
	w = 3;
	mode = 'r';

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
		if( strcmp( *(argv + i), "-n_lvl") == 0) {
			n_lvl  = atoi( *(argv + i + 1) );
			i++;
			free( nb_lvl );
			nb_lvl = (int *) malloc(n_lvl * sizeof(int));
 			for(j = 0; j < n_lvl; j++, i++) nb_lvl[j] = atoi( *(argv + i + 1) );
		}
		if( strcmp( *(argv + i), "-mt") == 0) {
			mt  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-mode") == 0) {
			if( strcmp( *(argv + i + 1), "levelx") == 0)    mode = 'l';
			if( strcmp( *(argv + i + 1), "recursive") == 0) mode = 'r';
			i++;
		}
		if( strcmp( *(argv + i), "-w") == 0) {
			w  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	if ( ( mode == 'l' )&&( w == 2 ) ) printf("dgeqrf_levelx_w02    | ");
	if ( ( mode == 'l' )&&( w == 3 ) ) printf("dgeqrf_levelx_w03    | ");
	if ( ( mode == 'r' )&&( w == 2 ) ) printf("dgeqrf_recursive_w02 | ");
	if ( ( mode == 'r' )&&( w == 3 ) ) printf("dgeqrf_recursive_w03 | ");

	printf("m = %4d, ",m);
	printf("n = %4d, ",n);
	printf("lda = %4d, ",lda);
	printf("ldq = %4d, ",ldq);
	if ( w == 3 ) printf("mt = %4d, ",mt); else printf("           "); 
//	printf("\n");
	if ( mode == 'l' ) {
	printf("n_lvl = %4d ( ",n_lvl);
 	for(j = 0; j < n_lvl; j++) printf(" %4d ",nb_lvl[j]);
	printf(")");
	if ( n_lvl == 1 ) {
	printf("            ");
	}
	if ( n_lvl == 2 ) {
	printf("      ");
	}
	}
	if ( mode == 'r' ) {
	printf("                                  ");
	}
	printf("  ");
//	printf("\n");

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

	if ( ( mode == 'l' )&&( w == 2 ) ){

	ldt = n;
	T = (double *) malloc(ldt * n * sizeof(double));

	int lwork;
	lwork = 1920000;
	work = (double *) malloc( 1920000 * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lila_dgeqrf_levelx_w02( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );
	free( T );
	}

	if ( ( mode == 'l' )&&( w == 3 ) ){

	ldt = mt;
	T = (double *) malloc(ldt * n * sizeof(double));

	int lwork;
	lwork = 1920000;
	work = (double *) malloc( 1920000 * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lila_dgeqrf_levelx_w03( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);


//	int k, i1;
//	for( k = 0; k < n; k++){
//	for( i1 = 0; i1 < n; i1++){
//	   printf(" %+5.2f ", T[ k + i1*ldt ] );
//	} 
//	printf("\n");
//	} 

	free( work );
	free( T );
	}

	if ( ( mode == 'r' )&&( w == 2 ) ){

	ldt = n;
	T = (double *) malloc(ldt * n * sizeof(double));

	int lwork;
	lwork = 1920000;
	work = (double *) malloc( 1920000 * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lila_dgeqrf_recursive_w02( m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );
	free( T );
	}

	if ( ( mode == 'r' )&&( w == 3 ) ){

	ldt = mt;
	T = (double *) malloc(ldt * n * sizeof(double));

	int lwork;
	lwork = 1920000;
	work = (double *) malloc( 1920000 * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lila_dgeqrf_recursive_w03( m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );
	free( T );
	}

	work = (double *) malloc(n * n * sizeof(double));

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );

	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );

	norm_orth = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );

	free( work );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, Q, ldq );

 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) As[i+j*lda] -= Q[i+j*ldq];

	normR = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, As, lda, work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;

	printf("| time = %f   GFlop/sec = %f   res = %e    orth = %e \n", elapsed_refL, perform_refL, normR / normA, norm_orth );

	free( nb_lvl );
	free( Q );
	free( A );
	free( As );
	free( S );

	return 0;
}
