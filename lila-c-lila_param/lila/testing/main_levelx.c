#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, ii, ml, verbose, testing;
	int lwork, n_lvl, *nb_lvl, panel, leaf, vrtq, *lila_param;
	double *A, *Q, *As, *T, *work=NULL, *Aii;
	double elapsed_ref1, perform_ref1, elapsed_ref2, perform_ref2;
	struct timeval tp;
	srand(0);

    	m         = 87;
    	n         = 53;
	ii        = 6;
	lda       = -1;
	ldq       = -1;
	mt        = 4;
	leaf      = 1;
	panel     = 1;
	n_lvl     = 1;
	nb_lvl    = (int *) malloc(n_lvl * sizeof(int));
	nb_lvl[0] = 10;
	verbose   = 0;
	testing   = 1;
	vrtq      = 0;


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
		if( strcmp( *(argv + i), "-ii") == 0) {
			ii  = atoi( *(argv + i + 1) );
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
	}

	if( m < n+ii ){ printf("\n\n YOUR CHOICE OF n AND ii HAVE MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	    	      ldt = mt;

	int k;
	k = 6;
	if( n_lvl > 1 ) k += n_lvl;
	lila_param = (int *) malloc(k * sizeof(int));
	lila_param[ 0 ] = 0;
	lila_param[ 1 ] = leaf;
	lila_param[ 2 ] = panel;
	lila_param[ 3 ] = n_lvl;
	lila_param[ 4 ] = vrtq;
	lila_param[ 5 ] = 0; // i_lvl starts at 0
	for(j = 0; j < n_lvl; j++) lila_param[6+j] = nb_lvl[j];

           	         A  = (double *) malloc(lda * (n+ii) * sizeof(double));
	                 As = (double *) malloc(lda * (n+ii) * sizeof(double));
	if ( vrtq != 1 ) Q  = (double *) malloc(ldq * (n+ii) * sizeof(double)); else Q = A;// if statement for wanting v03 don't allocate space for Q
   	                 T  = (double *) malloc(ldt * (n+ii) * sizeof(double));

 	for(i = 0; i < lda * (n+ii); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * (n+ii); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	ml   = m - ii;
	Aii  = A + ii + ii*lda;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+ii, A, lda, As, lda );

	lwork = lila_query_dgeqrf_w03_levelx( lila_param, n_lvl, 0, nb_lvl, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work=NULL, -1 );
	work  = (double *) malloc( lwork * sizeof(double));

	if( vrtq == 2 ){ 

		gettimeofday(&tp, NULL);
		elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_param[4] = 1; 
		lila_dgeqrf_levelx( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork ); 

		gettimeofday(&tp, NULL);
		elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		perform_ref1 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;

		gettimeofday(&tp, NULL);
		elapsed_ref2=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_param[4] = 2;
		lila_dgeqrf_levelx( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_ref2+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		perform_ref2 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref2 / 1.0e+9 ;

	} else if( vrtq == 3 ) { 

		gettimeofday(&tp, NULL);
		elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_param[4] = 1; 
		lila_dgeqrf_levelx( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork ); 

		gettimeofday(&tp, NULL);
		elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		perform_ref1 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;

		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', mt, n+ii, (0e+00), (0e+00), T, ldt );

		gettimeofday(&tp, NULL);
		elapsed_ref2=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_param[4] = 3;
		lila_dgeqrf_levelx( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_ref2+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		perform_ref2 = ( 0.0e+00 * ((double) m) * ((double) n) * ((double) n) - 0.0e+00 / 0.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref2 / 1.0e+9 ;

	} else {

		gettimeofday(&tp, NULL);
		elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_dgeqrf_levelx( lila_param, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
		perform_ref1 = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;

	}

	free( work );


	if ( verbose ){ 

		printf("dgeqrf_levelx    - ");
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
		printf("n_lvl = %4d ( ",n_lvl);
 		for(j = 0; j < n_lvl; j++) printf(" %4d ",nb_lvl[j]);
		printf(")");
		if ( n_lvl == 1 ) {
		printf("   ");
		}
		if ( n_lvl == 2 ) {
		printf("      ");
		}
		printf(" \n");

		printf(" time = %f    GFlop/sec = %f ", elapsed_ref1, perform_ref1);	
		if(( vrtq == 2 )||( vrtq == 3 )){ printf("    time2 = %f      Gflop/sec2 = %f ", elapsed_ref2, perform_ref2); }
		printf(" \n ");

	} else {

			printf("%6d %6d %6d %6d %6d %6d %6d %6d %6d ", m, ii, n, lda, ldq, mt, panel, leaf, n_lvl );
 			for(j = 0; j < n_lvl; j++) printf(" %4d ",nb_lvl[j]);
			printf(" %16.8f %10.3f ", elapsed_ref1, perform_ref1);
			if(( vrtq == 2 )||( vrtq == 3 )){ printf("%16.8f %10.3f ", elapsed_ref2, perform_ref2); }

	}

	if ( testing ){


		double orth_1, repres_1, repres_2, repres_3, repres_4;
	
//		vrtq == 0 || 2    ===>   we have Q to check
//		vrtq == 1 || 3    ===>   we don't
	
//		|| I - Q^T Q ||
		if (( vrtq == 0 )||( vrtq == 2 )) orth_1   = lila_test_qq_orth_1( m, n, ii, Q, ldq );

//		|| A - Q R || / || A ||
		if (( vrtq == 0 )||( vrtq == 2 )) repres_1 = lila_test_qr_repres_1( m, n, ii, As, lda, Q, ldq, A, lda );

//		|| (apply H) A - [R;0] || / || A ||
		repres_2 = lila_test_r_repres_2( lila_param, m, n, ii, mt, As, lda, T, ldt, A, lda );

//		create m-by-m H, then || I - HH^T ||; || H A - [R;0] || / || A ||
		repres_3 = lila_test_hh_repres( lila_param, m, n, ii, mt, As, lda, T, ldt, A, lda );

//		create m-by-n Q, then || I - Q^T Q ||; || A - Q R || / || A ||
		if (( vrtq == 1 )||( vrtq == 3 )) repres_4 = lila_test_vt_repres( lila_param, m, n, ii, mt, As, lda, T, ldt, A, lda );

		if (( vrtq == 0 )||( vrtq == 2 )){ if ( verbose ) printf("qq_orth   = %5.1e  \n ",orth_1); else printf(" %5.1e  ",orth_1); }
		if (( vrtq == 0 )||( vrtq == 2 )){ if ( verbose ) printf("qr_repres = %5.1e  \n ",repres_1); else printf(" %5.1e  ",repres_1); }
		if ( verbose ) printf("r_repres  = %5.1e  \n ",repres_2); else printf(" %5.1e  ",repres_2); 
		if ( verbose ) printf("h_q_orth  = %5.1e  \n ",repres_3); else printf(" %5.1e  ",repres_3); 
		if (( vrtq == 1 )||( vrtq == 3 )){ if ( verbose ) printf("vt_repres = %5.1e  \n ",repres_4); else printf(" %5.1e  ",repres_4); }

	}

	if ( !verbose ) printf("\n");		

	if ( vrtq != 1 ) free( Q );
	free( A );
	free( As );
	free( T );
	free( lila_param );

	return 0;
}
