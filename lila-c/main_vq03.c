#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, ii, lwork, ml, nx;
	int n_lvl, *nb_lvl, panel, leaf;
	int *lila_param;
	double *A, *Q, *As, *T, *work=NULL, *Aii;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;
	char mode;

	srand(0);

    	m         = 87;
    	n         = 53;
	ii        = 6;
	lda       = -1;
	ldq       = -1;
	mt        = 4;
	nx        = 7;
	leaf      = 0;
	panel     = 0;
	n_lvl     = 1;
	nb_lvl    = (int *) malloc(n_lvl * sizeof(int));
	nb_lvl[0] = 10;
	mode      = 'r';

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq  = atoi( *(argv + i + 1) );
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
		if( strcmp( *(argv + i), "-mt") == 0) {
			mt  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ii") == 0) {
			ii  = atoi( *(argv + i + 1) );
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
		if( strcmp( *(argv + i), "-n_lvl") == 0) {
			n_lvl  = atoi( *(argv + i + 1) );
			i++;
			free( nb_lvl );
			nb_lvl = (int *) malloc(n_lvl * sizeof(int));
 			for(j = 0; j < n_lvl; j++, i++) nb_lvl[j] = atoi( *(argv + i + 1) );
		}
		if( strcmp( *(argv + i), "-mode") == 0) {
			if( strcmp( *(argv + i + 1), "levelx") == 0)    mode = 'l';
			if( strcmp( *(argv + i + 1), "recursive") == 0) mode = 'r';
			i++;
		}
	}

	if( m < n+ii ){ printf("\n\n YOUR CHOICE OF n AND ii HAVE MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	if ( mode == 'l' ) printf("dgeqrf_levelx_w03    | ");
	if ( mode == 'r' ) printf("dgeqrf_recursive_w03 | ");

	printf("m = %4d, ",         m);
	printf("ii = %4d, ",       ii);
	printf("n = %4d, ",         n);
	printf("lda = %4d, ",     lda);
	printf("mt = %4d, ",       mt); 
	printf("panel = %4d, ", panel); 
	printf("leaf = %4d, ",   leaf); 
	if ( mode == 'l' ) {
	printf("n_lvl = %4d ( ",n_lvl);
 	for(j = 0; j < n_lvl; j++) printf(" %4d ",nb_lvl[j]);
	printf(")");
	if ( n_lvl == 1 ) {
	printf("   ");
	}
	if ( n_lvl == 2 ) {
	printf("      ");
	}
	}
	if ( mode == 'r' ) {
	printf("nx = %4d, ",       nx); 
	printf("                        ");
	}
	printf("  \n");

	A  = (double *) malloc(lda * (n+ii) * sizeof(double));
	As = (double *) malloc(lda * (n+ii) * sizeof(double));
	Q  = (double *) malloc(ldq * (n+ii) * sizeof(double));

 	for(i = 0; i < lda * (n+ii); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < lda * (n+ii); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	Aii   = A + ii + ii*lda;
	ml    = m - ii;
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A',  m, n+ii,   A, lda, As, lda );
	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml,    n, Aii, lda,    work );

	if ( mode == 'l' ){

		int k;
		k = 4;
		if( n_lvl > 1 ) k += n_lvl;
		lila_param = (int *) malloc(k * sizeof(int));
		lila_param[ 0 ] = 0;
		lila_param[ 1 ] = leaf;
		lila_param[ 2 ] = panel;
		lila_param[ 3 ] = n_lvl;

		if( lila_param[1] == 2 ){ printf("THE LU CODE DOES NOT WORK FOR v03"); return 0;}
		if( lila_param[2] == 2 ){ printf("THE LU CODE DOES NOT WORK FOR v03"); return 0;}

		ldt = mt;
		T = (double *) malloc(ldt * (n+ii) * sizeof(double));

		int lwork;
		work = NULL;
		lwork = lila_query_dgeqrf_w03_levelx( lila_param, n_lvl, 0, nb_lvl, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		work = (double *) malloc( lwork * sizeof(double));

		gettimeofday(&tp, NULL);
		elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_dgeqrf_v03_levelx( lila_param, n_lvl, 0, nb_lvl, m, n, ii, mt, A, lda, T, ldt, work, lwork );
		lila_dgeqrf_q03_levelx( lila_param, n_lvl, 0, nb_lvl, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		free( work );
	}

	if ( mode == 'r' ){

		lila_param = (int *) malloc(4 * sizeof(int));
		lila_param[ 0 ] = 0;
		lila_param[ 1 ] = leaf;
		lila_param[ 2 ] = panel;
		lila_param[ 3 ] = nx;

		if( lila_param[1] == 2 ){ printf("THE LU CODE DOES NOT WORK FOR v03"); return 0;}
		if( lila_param[2] == 2 ){ printf("THE LU CODE DOES NOT WORK FOR v03"); return 0;}

		ldt = mt;
		T = (double *) malloc(ldt * (n+ii) * sizeof(double));

		int lwork;
		work = NULL;
		lwork = lila_query_dgeqrf_w03_recursive( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		work = (double *) malloc( lwork * sizeof(double));

		gettimeofday(&tp, NULL);
		elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_dgeqrf_v03_recursive( lila_param, m, n, ii, mt, A, lda, T, ldt, work, lwork  );
		lila_dgeqrf_q03_recursive( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		free( work );
	}

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	
	info = lila_main_test( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, As, normA, elapsed_refL, perform_refL );

	free( A );
	free( Q );
	free( As );
	free( T );
	free( lila_param );

	return 0;
}
