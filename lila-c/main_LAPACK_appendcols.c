#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, vb, nb, verbose, testing, lwork;
	double *A, *Q, *As, *T, *work=NULL;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0);

    	m         = 87;
    	n         = 53;
	nb        = 10;
	lda       = -1;
	ldq       = -1;
	mt        = 1;
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
		if( strcmp( *(argv + i), "-verbose") == 0) {
			verbose  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-testing") == 0) {
			testing  = atoi( *(argv + i + 1) );
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
	}

	if( m < n ){ printf("\n\n YOUR CHOICE OF n HAS MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
		      ldt = 1;

	if ( verbose == 1 ){

		printf("dgeqrf_w03_recursive | ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("lda = %4d, ",     lda);
		printf("ldq = %4d, ",     ldq);
		printf("mt = %4d, ",       mt);
		}

	A  = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q  = (double *) malloc(ldq * n * sizeof(double));
	T  = (double *) malloc( 1  * n * sizeof(double));

 	for(i = 0; i < lda * (n); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * (n); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	normA = LAPACKE_dlange_work   ( LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

	j = 0;
	if ( nb > n ) vb = n; else vb = nb;
	lwork = 0;
	while( vb!=0 ){

		lwork = lila_query_dgeqrf_LAPACK_appendcols( m, j, vb, 1, A, lda, T, ldt, Q, ldq, work=NULL, lwork );
		j += vb;
		if ( j+nb > n ) vb = n-j; else vb = nb;

	}

	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	j = 0;
	if ( nb > n ) vb = n; else vb = nb;
	while( vb!=0 ){

		info = lila_dgeqrf_LAPACK_appendcols( m, j, vb, 1, A, lda, Q, ldq, T, ldt, work, lwork );
		j += vb;
		if ( j+nb > n ) vb = n-j; else vb = nb;

	}

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	
	if ( verbose == 0 ){ 
		printf("%6d %6d %6d  %16.8f %10.3f \n", m, n, mt, elapsed_refL, perform_refL);
		//printf("%6d %6d %6d  %s %16.8f %10.3f \n", m, n, mt, getenv("OPENBLAS_NUM_THREADS"), elapsed_refL, perform_refL);
	} 

	if ( testing == 1 ){
		int *lila_param;
		lila_param = (int *) malloc(4 * sizeof(int));
		lila_param[ 0 ] = 0;
		lila_param[ 1 ] = 0;
		lila_param[ 2 ] = 0;

		info = lila_main_test( lila_param, m, n, 0, mt, A, lda, T, ldt, Q, ldq, As, normA, elapsed_refL, perform_refL );

		free( lila_param );
	}

	free( Q );
	free( A );
	free( As );
	free( T );

	return 0;
}
