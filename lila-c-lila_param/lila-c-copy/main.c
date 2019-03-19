#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, ii, ml, nx, verbose, testing;
	int lwork, n_lvl, *nb_lvl, panel, leaf, vrtq, t03, *lila_param;
	double *A, *Q, *As, *T, *Ts, *work=NULL, *Aii;
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
	leaf      = 1;
	panel     = 1;
	n_lvl     = 1;
	nb_lvl    = (int *) malloc(n_lvl * sizeof(int));
	nb_lvl[0] = 10;
	mode      = 'r';
	verbose   = 0;
	testing   = 0;
	vrtq      = 0;
	t03	  = 1;


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
		if( strcmp( *(argv + i), "-mt") == 0) {
			mt  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-vrtq") == 0) {
			vrtq  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-t03") == 0) {
			t03  = atoi( *(argv + i + 1) );
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
	    	      ldt = mt;

	if ( verbose == 1 ){

		if ( mode == 'r' ) printf("dgeqrf_recursive - ");
		if ( mode == 'l' ) printf("dgeqrf_levelx    - ");
		if ( vrtq == 0 ) printf(" w03 | ");
		if ( vrtq == 1 ) printf(" v03 | ");
		if ( vrtq == 2 ) printf(" v03 & q03 | ");
		if ( vrtq == 3 ) printf(" t03 | ");
		printf("m = %4d, ",         m);
		printf("ii = %4d, ",       ii);
		printf("n = %4d, ",         n);
		printf("lda = %4d, ",     lda);
		printf("ldq = %4d, ",     ldq);
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
		printf("nx = %4d, ", nx); 
		printf("                        ");
		}
		printf("  ");
	}

	if ( mode == 'r' ){
		if ( vrtq == 3 ){
			lila_param = (int *) malloc(6 * sizeof(int));
			lila_param[ 0 ] = mode;
			lila_param[ 1 ] = leaf;
			lila_param[ 2 ] = panel;
			lila_param[ 3 ] = nx;
			lila_param[ 4 ] = vrtq;
			lila_param[ 5 ] = t03;
		} else {
			lila_param = (int *) malloc(5 * sizeof(int));
			lila_param[ 0 ] = mode;
			lila_param[ 1 ] = leaf;
			lila_param[ 2 ] = panel;
			lila_param[ 3 ] = nx;
			lila_param[ 4 ] = vrtq;
		}
	} else {
		int k;
		k = 5;
		if( n_lvl > 1 ) k += n_lvl;
		lila_param = (int *) malloc(k * sizeof(int));
		lila_param[ 0 ] = mode;
		lila_param[ 1 ] = leaf;
		lila_param[ 2 ] = panel;
		lila_param[ 3 ] = n_lvl;
		lila_param[ 4 ] = vrtq;
		for(j = 0; j < n_lvl; j++) nb_lvl[j] = lila_param[5+j];

	}


           	         A  = (double *) malloc(lda * (n+ii) * sizeof(double));
	                 As = (double *) malloc(lda * (n+ii) * sizeof(double));
	if ( vrtq != 1 ) Q  = (double *) malloc(ldq * (n+ii) * sizeof(double)); else Q = A;// if statement for wanting v03 don't allocate space for Q
   	                 T  = (double *) malloc(ldt * (n+ii) * sizeof(double));
	if ( vrtq == 3 ) Ts = (double *) malloc( n  *  n     * sizeof(double));

 	for(i = 0; i < lda * (n+ii); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * (n+ii); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	ml   = m - ii;
	Aii  = A + ii + ii*lda;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+ii, A, lda, As, lda );

	lwork = lila_query_dgeqrf_w03_recursive( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work=NULL, -1 );
	work  = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	////
	if( vrtq == 2 ){ lila_param[4] = 1; lila_dgeqrf( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork ); lila_param[4] = vrtq;}
	lila_dgeqrf( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	////

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	
	if ( verbose == 0 ){ 
		if ( mode == 'r' ){
			printf("%6d %6d %6d %6d %16.8f %10.3f %6d %6d\n", m, n, mt, nx, elapsed_refL, perform_refL, leaf, panel);
			//printf("%6d %6d %6d %6d %s %16.8f %10.3f %6d %6d\n", m, n, mt, nx, getenv("OPENBLAS_NUM_THREADS"), elapsed_refL, perform_refL, leaf, panel);
		} else { // levelx
			printf("%6d %6d %6d %16.8f %10.3f %6d %6d\n", m, n, mt, elapsed_refL, perform_refL, leaf, panel);
			//printf("%6d %6d %6d %s %16.8f %10.3f %6d %6d\n", m, n, mt, getenv("OPENBLAS_NUM_THREADS"), elapsed_refL, perform_refL, leaf, panel);
		} 
	} 

	if ( testing == 1 ){
		normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, n, As+ii*(1+lda), lda, work );
		if (( vrtq == 1 ) || ( vrtq == 3 )){
			info = lila_main_test( lila_param, m, n, ii, mt, A, lda, T, ldt, NULL, -1, As, normA, elapsed_refL, perform_refL );
		} else {
			info = lila_main_test( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, As, normA, elapsed_refL, perform_refL );
		}
/*		
		// This if statement compares the T constructed using our dlarft_w03, LAPACKE_dlarft_work (with the mt-structure), with the nxn Ts  
		// After t03 is working we can move this within its own script
		if ( vrtq == 3 ){

			double norm_diff_T, *Tjj;
			int vb, jj;
			double *tau, *Akk, normv2;
			int k; 

			tau = (double *) malloc(n * sizeof(double));
			Akk=Aii;
			info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (0e+00), Ts, n );// this line is not needed. Just being safe
			for( k = 0; k < n; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }
			info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Ts, n);
	
			//for( k = 0; k < n; k++ ) printf("T = %2.3e, Ts = %2.3e \n", T[(k%mt)+k*ldt], Ts[k+k*n]); // print the diagonal (tau)	

			jj = ii;
			Tjj = T+jj+jj*ldt;
			vb = mt - ( ii%mt ); if ( vb > n ) vb = n;

	 		for(i = 0; i < vb; i++) for(j = 0; j < vb; j++) Tjj[ i+j*ldt ] -= Ts[ i+j*n ];			

			jj += vb;
			Tjj = T+jj+jj*ldt;
			if( jj + mt >= ii + n ) vb = n - ( jj - ii ); else vb = mt;

			while( vb != 0 ){

		 		for(i = 0; i < vb; i++) for(j = jj; j < jj+vb; j++) Tjj[ i+j*ldt ] -= Ts[ (i+jj)+j*n ];

				jj += vb;
				Tjj = T+jj+jj*ldt;
				if( jj + mt >= ii + n ) vb = n - ( jj - ii ); else vb = mt;

			}

			norm_diff_T = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', mt, n, T, ldt, NULL );
			printf("  diff_t = %2.1e \n", norm_diff_T );
			free( tau ) ;
		}
*/
	}


	if ( vrtq != 1 ) free( Q );
	free( A );
	free( As );
	if ( vrtq == 3 ) free( Ts );
	free( T );
	free( lila_param );

	return 0;
}
