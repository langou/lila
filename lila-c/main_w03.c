#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt, ii, lwork, ml, nx;
	int n_lvl, *nb_lvl, panel, leaf;
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
	leaf      = 1;
	panel     = 1;
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

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	if ( mode == 'l' ) printf("dgeqrf_levelx_w03    | ");
	if ( mode == 'r' ) printf("dgeqrf_recursive_w03 | ");

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
	printf("nx = %4d, ",       nx); 
	printf("                        ");
	}
	printf("  ");

	A  = (double *) malloc(lda * (n+ii) * sizeof(double));
	As = (double *) malloc(lda * (n+ii) * sizeof(double));
	Q  = (double *) malloc(ldq * (n+ii) * sizeof(double));

 	for(i = 0; i < lda * (n+ii); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * (n+ii); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	Aii   = A + ii + ii*lda;
	ml    = m-ii;
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+ii, A, lda, As, lda );
	normA = LAPACKE_dlange_work   ( LAPACK_COL_MAJOR, 'F', ml, n, Aii, lda, work );

	if( m < n+ii ){ printf("\n\n YOUR CHOICE OF n AND ii HAVE MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if ( mode == 'l' ){

		ldt = mt;
		T = (double *) malloc(ldt * (n+ii) * sizeof(double));

		int lwork;
		lwork = 0;
		work = (double *) malloc( 1920000 * sizeof(double));

		lwork = lila_wsq_dgeqrf_levelx_w03( panel, leaf, n_lvl, 0, nb_lvl, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		//printf(" 0 |  lwork  = %3d,\n",lwork);
		free( work );
		work = (double *) malloc( lwork * sizeof(double));
		//work = (double *) malloc( 1920000 * sizeof(double));

		gettimeofday(&tp, NULL);
		elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_dgeqrf_levelx_w03( panel, leaf, n_lvl, 0, nb_lvl, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		free( work );
	}

	if ( mode == 'r' ){

		ldt = mt;
		T = (double *) malloc(ldt * (n+ii) * sizeof(double));

		int lwork;
		lwork = 0;
		work = (double *) malloc( 1920000 * sizeof(double));

		lwork = lila_wsq_dgeqrf_recursive_w03( panel, leaf, nx, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		//printf(" 0 |  lwork  = %3d,\n",lwork);
		free( work );
		work = (double *) malloc( lwork * sizeof(double));
		//work = (double *) malloc( 1920000 * sizeof(double));

		gettimeofday(&tp, NULL);
		elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		lila_dgeqrf_recursive_w03( panel, leaf, nx, m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

		gettimeofday(&tp, NULL);
		elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

		free( work );
	}

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	
	double *QQ, *RR, *HH, norm_repres_2_1, norm_repres_2_2, norm_orth_2, *Qii, *Tii, *As_ii;
	double norm_orth_3, norm_repres_3, norm_diffQ_3, norm_orth_1, norm_repres_1;

	As_ii = As+ii+ii*lda;
	Qii   =  Q+ii+ii*ldq;
	Tii   =  T+ii+ii*ldt;

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
	if ( ( leaf == 1 )&&( panel == 2) ){
		lila_dormqrf_z03( m, n, n, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	} else if( leaf == 2 ) {
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
	if( ( leaf == 1 )&&( panel == 2) ){
		lila_dormqrf_z00( ml, ml, n,  0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	} else if ( leaf == 2 ) {
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


	free( Q );
	free( A );
	free( As );
	free( T );

	return 0;
}
