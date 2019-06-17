#include "ls.h"

int main(int argc, char ** argv) {

	int i, m, n, k, lda, ldr, ldb, ldx, lwork, info, verbose, testing;
	double *A, *As, *R, *x, *b, *work;
	double norm_A, norm_b; 
	double norm_repres, norm_orth;
	double elapsed, perform;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 21;
    	k         =  1;
	lda       = -1;
	ldr       = -1;
	ldb       = -1;
	ldx       = -1;
	verbose   = 0;
	testing   = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldr") == 0) {
			ldr = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldb") == 0) {
			ldb = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldx") == 0) {
			ldx = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-verbose") == 0) {
			verbose = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-testing") == 0) {
			testing = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-m") == 0) {
			m = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-k") == 0) {
			k = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if (( m < n )) { printf("\n\n We need n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldr < 0 ) ldr = n;
	if( ldb < 0 ) ldb = m;
	if( ldx < 0 ) ldx = m;

	A  = (double *) malloc( lda * n * sizeof(double));
	As = (double *) malloc( lda * n * sizeof(double));
 	R  = (double *) malloc( ldr * n * sizeof(double));
	x  = (double *) malloc( ldx * k * sizeof(double));
	b  = (double *) malloc( ldb * k * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldr * n; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldx * k; i++)
		*(x + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, x, ldx,  b, ldb );
	norm_A = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, As, lda, NULL );
	norm_b = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k,  b, ldb, NULL );
	

	work = (double *) malloc( 1 * sizeof(double));
	LAPACKE_dgels_work( LAPACK_COL_MAJOR, 'N', m, n, k, A, lda, x, ldx, work, -1 ); 
	lwork = ((int) work[0]);
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	info = LAPACKE_dgels_work( LAPACK_COL_MAJOR, 'N', m, n, k, A, lda, x, ldx, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	perform = ( 1.00 ) / elapsed / 1.0e+9 ; // need the flop count for gels -- flop count of geqrf and orgqr + other things

	if ( verbose ){ 

		printf("LAPACK dGELS");
		printf("m = %4d, ", m);
		printf("n = %4d, ", n);
		printf("k = %4d, ", k);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed, perform );	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %16.8f %10.3f ", m, n, k, elapsed, perform );

	} 

	if ( testing ){


		int ii, jj;
		double *r;
		double norm_lower_x, norm_x, *orth_vec;
		orth_vec = (double *) malloc( n * k * sizeof(double));
		r        = (double *) malloc( m * k * sizeof(double));

		norm_x = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, k, x, ldx, NULL );	
		norm_lower_x = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m-n, k, x+n, ldx, NULL );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, n, (+1.0e+00), As, lda, x, ldx, (+0.0e+00), r, m );
 		for(jj = 0; jj < k; jj++) for(ii = 0; ii < m; ii++) r[ ii+jj*m ] -= b[ ii+jj*ldb ];
		norm_repres = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, r, m, NULL );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n, k, m, (+1.0e+00), As, lda, r, m, (+0.0e+00), orth_vec, n );
		norm_orth = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, k, orth_vec, n, NULL );

		if ( verbose ) printf("norm_orth = %3.3e ", norm_orth / ( norm_A * ( norm_A * norm_x + norm_b ) ) );
		else printf( "%3.3e ", norm_orth / ( norm_A * ( norm_A * norm_x + norm_b ) ) );
		if ( verbose ) printf("norm_repres = %3.3e ", ( norm_repres -  norm_lower_x ) / ( norm_A * norm_x + norm_b ) );
		else printf(" %3.3e ", ( norm_repres -  norm_lower_x ) / ( norm_A * norm_x + norm_b ) );

		free( r );
		free( orth_vec );

	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( As );
	free( R  );
	free( b  );
	free( work );

	return 0;

}
