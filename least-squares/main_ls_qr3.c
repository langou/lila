#include "ls.h"

int main(int argc, char ** argv) {

	int i, j, m, n, k, lda, ldt, ldb, verbose, testing;
	double *A, *As, *T, *bs, *b;
	double norm_A, norm_b; 
	double norm_repres, norm_orth;
	double elapsed, perform;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 21;
    	k         =  1;
	lda       = -1;
	ldt       = -1;
	ldb       = -1;
	verbose   = 0;
	testing   = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldt") == 0) {
			ldt = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldb") == 0) {
			ldb = atoi( *(argv + i + 1) );
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

	if (( m < n ) && ( n < k )) { printf("\n\n We need k <= n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldb < 0 ) ldb = m;
	if( ldt < 0 ) ldt = n;

	A  = (double *) malloc( lda * n * sizeof(double));
	T  = (double *) malloc( ldt * n * sizeof(double));
	As = (double *) malloc( lda * n * sizeof(double));
	bs = (double *) malloc( ldb * k * sizeof(double));
	b  = (double *) malloc( ldb * k * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldb * k; i++)
		*(b + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, b, ldb, bs, ldb );
	norm_A = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, As, lda, NULL );
	norm_b = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, bs, ldb, NULL );
	
//	double *work;
//	int lwork;
//	work = (double *) malloc( 1 * sizeof(double));
//	LAPACKE_dgels_work( LAPACK_COL_MAJOR, 'N', m, n, k, A, lda, x, ldx, work, -1 ); 
//	lwork = ((int) work[0]);
//	free( work );
//	work = (double *) malloc( lwork * sizeof(double));
//	LAPACKE_dgels_work( LAPACK_COL_MAJOR, 'N', m, n, k, A, lda, x, ldx, work, lwork );
//	free( work );

	gettimeofday(&tp, NULL);
	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	dgeqr3_R( m, n, A, lda, T, ldt, A, lda );

	double *work;
	work = (double *) malloc( n * k * sizeof(double));

	// copy b into work
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, k, b, ldb, work, n ); 

	// work = V1^T * b1
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n, k, (+1.0e+00), A, lda, work, n );
 
	// work = V2^T * b2 + work
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n, k, m-n, (+1.0e+00), A+n, lda, b+n, ldb, (+1.0e+00), work, n );

	// work = T * work
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n, k, (+1.0e+00), T, ldt, work, n ); 

	// x2 = b2 - V2 * work 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n, k, n, (+1.0e+00), A+n, lda, work, n, (+0.0e+00), b+n, ldb );
	for(j=0;j<k;j++) for(i=0;i<m-n;i++) b[ i+n+j*ldb ] = bs[ i+n+j*ldb ] - b[ i+n+j*ldb ];

	// x1 = b1 - V1 * work
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, k, (+1.0e+00), A, lda, work, n ); 
	for(j=0;j<k;j++) for(i=0;i<n;i++) b[ i+j*ldb ] = bs[ i+j*ldb ] - work[ i+j*n ];

	// Solve for x
	cblas_dtrsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, k, (+1.0e00), A, lda, b, ldb );

	free( work );

	gettimeofday(&tp, NULL);
	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	perform = ( 
		2.0e+00 * ( (double) m) * ((double) n) *((double) n) - ((double) n) * ((double) n) * ((double) n) / 3.0e+00 
		+ 4.0e+00 * ((double) (m-n)) * ((double) n) * ((double) k)
		+ 4.0e+00 * ((double) n) * ((double) n) * ((double) k) 
		) / elapsed / 1.0e+9 ; // need the flop count for gels -- flop count of geqrf and orgqr + other things

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

		norm_x = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, k, b, ldb, NULL );	
		norm_lower_x = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m-n, k, b+n, ldb, NULL );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m, k, n, (+1.0e+00), As, lda, b, ldb, (+0.0e+00), r, m );
 		for(jj = 0; jj < k; jj++) for(ii = 0; ii < m; ii++) r[ ii+jj*m ] -= bs[ ii+jj*ldb ];
		norm_repres = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, r, m, NULL );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n, k, m, (+1.0e+00), As, lda, r, m, (+0.0e+00), orth_vec, n );
		norm_orth = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, k, orth_vec, n, NULL );

		if ( verbose ) printf("norm_orth = %3.3e ", norm_orth / ( norm_A * ( norm_A * norm_x + norm_b ) ) );
		else printf(" %3.3e ", norm_orth / ( norm_A * ( norm_A * norm_x + norm_b ) ) );
		if ( verbose ) printf("norm_repres = %3.3e ", ( norm_repres -  norm_lower_x ) / ( norm_A * norm_x + norm_b ) );
		else printf(" %3.3e ", ( norm_repres -  norm_lower_x ) / ( norm_A * norm_x + norm_b ) );

		free( r );
		free( orth_vec );

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( T );
	free( bs );
	free( b );
	free( As );

	
	return 0;

}
