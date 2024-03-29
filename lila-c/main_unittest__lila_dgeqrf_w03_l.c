#include "lila.h"

int main(int argc, char ** argv){

	int m, n, mt, nb, ml, ii, i, j, lda, ldq, ldt, info, lwork;
	double *A, *As, *T, *Q, *work=NULL;
	double normA, elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0);

    	m   = 41;
    	n   = 23;
	nb  =  2;
	ii  =  1;
	mt  =  6;
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
		if( strcmp( *(argv + i), "-ii") == 0) {
			ii  = atoi( *(argv + i + 1) );
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
	printf("ii = %4d; ",ii);
	printf("mt = %4d; ",mt);
	printf("nb = %4d; ",nb);
	printf("lda = %4d; ",lda);
	printf("ldq = %4d; ",ldq);
	printf("  ");

	A = (double *) malloc(lda * (n+ii) * sizeof(double));
	As = (double *) malloc(lda * (n+ii) * sizeof(double));
	Q = (double *) malloc(ldq * (n+ii) * sizeof(double));

 	for(i = 0; i < lda * (n+ii); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
 	for(i = 0; i < ldq * (n+ii); i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	double *Aii;
	Aii = A+ii+ii*lda;
	ml = m-ii;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+ii, A, lda, As, lda );
	normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', ml, n, Aii, lda, work );

	ldt = mt;
	T = (double *) malloc(ldt * (n+ii) * sizeof(double));

	printf("\n");
	lwork = 200000000;
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lila_dgeqrf_w03_level1( m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork, nb );

//	lila_dgeqrf_w03_mt    ( m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//	lila_dgeqrf_w03_l     ( m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

//	lila_dgeqrf_w03_mt_v02( m, n, ii, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

/*
	double *RRR;
	double *RRRi0, *AAs_ii;
	RRR = (double *) malloc(m * (n+ii) * sizeof(double));
	RRRi0 = RRR+ii;
	AAs_ii = As+ii+ii*lda;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, AAs_ii, lda, RRRi0, m );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, AAs_ii, lda, RRRi0, m );
	lila_dormqrf_z03( m, n, n, ii, 0, mt, A, lda, T, ldt, RRR, m, work, lwork );
	//info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, RRRi0, lda, Aii, lda ); 

	int ik;
	int jk;
	printf("R=[\n");
	for( ik = ii; ik < n+ii; ik++ ){
		for( jk = ii; jk < ik; jk++ ){
			printf("          ");
		}
		for( jk = ik; jk < n+ii; jk++ ){
			//printf(" %6.2f ", A[ ik + jk*lda]-RRR[ ik + (jk-ii)*lda]);
			if (fabs(A[ ik + jk*lda]-RRR[ ik + (jk-ii)*lda])>0.00001){
			printf(" *%6.2f* ", A[ ik + jk*lda]);
			}else{
			printf("  %6.2f  ", A[ ik + jk*lda]);
			}
		}
		printf("\n");
	}
	printf("];\n");

	printf("R=[\n");
	for( ik = ii; ik < n+ii; ik++ ){
		for( jk = 0; jk < ik-ii; jk++ ){
			printf("          ");
		}
		for( jk = ik-ii; jk < n; jk++ ){
			if (fabs(A[ ik + (jk+ii)*lda]-RRR[ ik + jk*lda])>0.00001){
			printf(" *%6.2f* ", RRR[ ik + jk*lda]);
			}else{
			printf("  %6.2f  ", RRR[ ik + jk*lda]);
			}
		}
		printf("\n");
	}
	printf("];\n");


	free( RRR );
*/

	free( work );

	perform_refL = ( 4.0e+00 * ((double) ml) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;

////////////////////////////////////////////////////////////////////////////////////////////////////
	double norm_orth_1, norm_repres_1;
	double *QQ, *RR, *HH, norm_repres_2_1, norm_repres_2_2, norm_orth_2;
	double norm_orth_3, norm_repres_3, norm_diffQ_3;
	double *Qii, *Tii, *As_ii;

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

//	RR = (double *) malloc(ml * n * sizeof(double));
	RR = (double *) malloc(m * n * sizeof(double));
	double *RRi0;
	RRi0 = RR+ii;
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, As_ii, lda, RRi0, m );

	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
//	if( ml > n )  lila_dormqrf_z03( m, n, n,   ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
//	if( ml == n ) lila_dormqrf_z03( m, n, n-1, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
	lila_dormqrf_z03( m, n, n, ii, 0, mt, A, lda, T, ldt, RR, m, work, lwork );
//
//	if( ml > n )  lila_dormqrf_z03( ml, n, n,   0, 0, mt, A+ii*(1+lda), lda, T+ii*(1+ldt), ldt, RR, ml, work, lwork );
//	if( ml == n ) lila_dormqrf_z03( ml, n, n-1, 0, 0, mt, A+ii*(1+lda), lda, T+ii*(1+ldt), ldt, RR, ml, work, lwork );
	free( work );
//	printf("3 \n");

//	if( ml > n )  norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', ml-1, n,   RRi0+1, m, NULL );
//	if( ml == n ) norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', ml-1, n-1, RRi0+1, m, NULL );
	norm_repres_2_1 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'L', 'N', ml-1, n, RRi0+1, m, NULL );

	lwork = n*n;
	work  = (double *) malloc(n * n * sizeof(double));
 //	for(i = 0; i < n; i++) for(j = 0; j < n; j++) work[ i+j*n ] = A[ (i+ii)+(j+ii)*lda] - RR[ i+j*ml ];
 	for(i = 0; i < n; i++) for(j = 0; j < n; j++) work[ i+j*n ] = Aii[ i+j*lda ] - RRi0[ i+j*m ];
	norm_repres_2_2 = LAPACKE_dlantr_work(LAPACK_COL_MAJOR, 'F', 'U', 'N', n, n, work, n, NULL );
	free( work );

	free( RR );
//	printf("4 \n");

	HH = (double *) malloc(ml * ml * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', ml, ml, (0e+00), (1e+00), HH, ml );
	lwork = ml*ml;
	work  = (double *) malloc(ml * ml * sizeof(double));
//	if( ml > n )  lila_dormqrf_z03( ml, ml, n,   0, 0, mt, Aii, lda, HH, ml, Tii, ldt, work, lwork );
//	if( ml == n ) lila_dormqrf_z03( ml, ml, n-1, 0, 0, mt, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	if( ml > n )  lila_dormqrf_z00( ml, ml, n,   0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
	if( ml == n ) lila_dormqrf_z00( ml, ml, n-1, 0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
//	lila_dormqrf_z00( ml, ml, n, 0, 0, -1, Aii, lda, HH, ml, Tii, ldt, work, lwork );
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
