#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, nx, vb, nb, verbose, testing;
	int panel, leaf, *lila_param, lwork;
	double *A, *Q, *As, *T, *work=NULL;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0);

    	m         = 87;
    	n         = 53;
	nb        = 10;
	lda       = -1;
	ldq       = -1;
	mt        = 4;
	nx        = 1;
	leaf      = 1;
	panel     = 1;
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
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-mt") == 0) {
			mt  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-nx") == 0) {
			nx  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-panel") == 0) {
			panel  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-leaf") == 0) {
			leaf  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( m < n ){ printf("\n\n YOUR CHOICE OF n HAS MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
		      ldt = mt;

	lila_param = (int *) malloc(4 * sizeof(int));
	lila_param[ 0 ] = 0;
	lila_param[ 1 ] = leaf;
	lila_param[ 2 ] = panel;
	lila_param[ 3 ] = nx;

	if ( verbose == 1 ){

		printf("dgeqrf_w03_recursive | ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("nb = %4d, ",       nb);
		printf("lda = %4d, ",     lda);
		printf("ldq = %4d, ",     ldq);
		printf("mt = %4d, ",       mt);
		printf("panel = %4d, ", panel); 
		printf("leaf = %4d, ",   leaf); 
		printf("nx = %4d, ",       nx); 
		}

	A  = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q  = (double *) malloc(ldq * n * sizeof(double));
	T  = (double *) malloc(ldt * n * sizeof(double));

 	for(i = 0; i < lda * (n); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * (n); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	normA = LAPACKE_dlange_work   ( LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

		lwork = lila_query_dgeqrf_w03_recursive( lila_param, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work=NULL, -1 );
		work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	j = 0;
	if ( nb > n ) vb = n; else vb = nb;
	while( vb!=0 ){

		lila_dgeqrf_w03_appendcols( lila_param, m, j, vb, mt, A, lda, Q, ldq, T, ldt, work, lwork );
		j += vb;
		if ( j+nb > n ) vb = n-j; else vb = nb;

	}

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	
	if ( verbose == 0 ){ 
		printf("%6d %6d %6d %6d %6d %16.8f %10.3f %6d %6d\n", m, n, nb, mt, nx, elapsed_refL, perform_refL, leaf, panel);
		//printf("%6d %6d %6d %6d %s %16.8f %10.3f %6d %6d\n", m, n, mt, nx, getenv("OPENBLAS_NUM_THREADS"), elapsed_refL, perform_refL, leaf, panel);
	} 

	if ( testing == 1 ){
		info = lila_main_test( lila_param, m, n, 0, mt, A, lda, T, ldt, Q, ldq, As, normA, elapsed_refL, perform_refL );
	}

	free( Q );
	free( A );
	free( As );
	free( T );
	free( lila_param );

	return 0;
}
