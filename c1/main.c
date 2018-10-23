#include "lila.h"

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, ldt, m, n, mt;
	double *A, *Q, *As, *T, *work=NULL;
	double normA, normR;
	double elapsed_refL, perform_refL;
	struct timeval tp;
	int n_lvl;
	int *nb_lvl;

	srand(0);

    	m = 30;
    	n = 17;
	lda = -1;
	ldq = -1;

	mt = 9;
	n_lvl = 1;
	nb_lvl = (int *) malloc(n_lvl * sizeof(int));
	nb_lvl[0] = 10;

	double *TTT;
	int llldddttt;

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

	}

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	ldt = n;

	printf("m = %4d, ",m);
	printf("n = %4d, ",n);
	printf("lda = %4d, ",lda);
	printf("ldq = %4d, ",ldq);
	printf("mt = %4d, ",mt);
	printf("\n");
	printf("n_lvl = %4d ( ",n_lvl);
 	for(j = 0; j < n_lvl; j++) printf(" %4d ",nb_lvl[j]);
	printf(")");
	printf("\n");

	A = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q = (double *) malloc(ldq * n * sizeof(double));
	T = (double *) malloc(ldt * n * sizeof(double));

//	llldddttt = ldt;	
	llldddttt = mt;	
	TTT = (double *) malloc(llldddttt * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//	lila_dge_qr_wq_vL0( m, n, 0, mt, A, lda, T, ldt, Q, ldq );
//	lila_dge_qr_wq_vr0( m, n, 0, mt, A, lda, T, ldt, Q, ldq );
//	lila_dge_qr_wq_manylevels_vr0( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq );

	int lwork;
	lwork = 1920000;
	work = (double *) malloc( 1920000 * sizeof(double));

//lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );

//lila_dge_qr_wq_levelx_w00( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//lila_dge_qr_wq_levelx_w02( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//lila_dge_qr_wq_levelx_w03( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );
lila_dge_qr_wq_levelx_w03( n_lvl, 0, nb_lvl, m, n, 0, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork );

//lila_dgeqr2_recursive( m, n, 0, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork );


free(TTT);

	free(work);

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, Q, ldq );

 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) As[i+j*lda] -= Q[i+j*ldq];

	normR = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, As, lda, work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;

	printf("LAPACK        :: time = %f   GFlop/sec = %f   res = %e \n", elapsed_refL, perform_refL, normR / normA );

	free( nb_lvl );
	free( T );
	free( Q );
	free( A );
	free( As );

	return 0;
}
