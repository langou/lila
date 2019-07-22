#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int check_qq_orth( double *norm_orth_1, int m, int n, double *Q, int ldq  );
extern int check_qr_repres( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );

extern int dlatsqr_( int *m, int *n, int *mb, int *nb, double *A, int *lda, double *T, int *ldt, double *work, int *lwork, int *info);

extern int dlamtsqr_( char *side, char *trans, int *m, int *n, int *k, int *mb, int *nb,
		double *A, int *lda, double *T, int *ldt, double *C, int *ldc, double *work, int *lwork, int *info );

//extern int dtpmqrt_( char *side, char *trans, int *m, int *n, int *k, int *l, int *nb, 
//		double *V, int *ldv, double *T, int *ldt, double *A, int *lda, double *B, int *ldb, double *work, int *info );

//extern int dgemqrt_( char *side, char *trans, int *m, int *n, int *k, int *nb,
//		double *v, int *ldv, double *t, int *ldt, double *c, int *ldc, double *work, int *info );

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, ldt, m, n, verbose, testing;
	int lwork;
	double *A, *Q, *R, *T, *tau, *work, *Asave;
	double orth, repres;
	double elapsed, perform;
	struct timeval tp;
	int mb, nb;
	
	srand(0);

    	m         = 81;
    	n         = 20;
	mb        = 25;
	nb        =  5;
	lda       = -1;
	ldq       = -1;
	ldr       = -1;
	ldt       = -1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldr") == 0) {
			ldr = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldt") == 0) {
			ldt = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-m") == 0) {
			m = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-mb") == 0) {
			mb = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if ( m < n ) { printf("\n\n We need n <= m\n\n"); return 0; }
	if ( mb < n ) { printf("\n\n We need n <= mb\n\n"); return 0; }
	if ( n < nb ) { printf("\n\n We need nb <= n\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	if( ldr < 0 ) ldr = n;
	if( ldt < 0 ) ldt = n;

	A = (double *) malloc( lda * n * sizeof(double));
	Asave = (double *) malloc( lda * n * sizeof(double));
	Q = (double *) malloc( ldq * n * sizeof(double));
 	R = (double *) malloc( ldr * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < lda * n; i++)
		*(Asave + i) = *(A + i);

	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldr * n; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

/*
	tau = (double *) malloc( n * sizeof(double));

	work = (double *) malloc( 1 * sizeof(double));
	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, -1 ); 
	lwork = ((int) work[0]);
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, -1 );
	if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );
	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, lwork ); 
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Q, ldq, R, ldr );
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( tau );
	free( work );
*/


/*
	int info;

 	T = (double *) malloc( 10 * ldt * n * ((int) ceil(( (double) (m-n) )/((double) (mb-n)))) * sizeof(double));

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );
	lwork = ((int) work[0]);
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, A, lda, R, ldr );

	LAPACKE_dlaset_work( LAPACK_COL_MAJOR, 'A', m, n, 0.0e+00, 1.0e+00, Q, ldq );

	dlamtsqr_( "L", "N", &m, &n, &n, &mb, &nb, A, &lda, T, &ldt, Q, &ldq, work, &lwork, &info );

	free( work );
	free( T );
*/



/*
	int info;

 	T = (double *) malloc( ldt * n * ((int) ceil(( (double) (m-n) )/((double) (mb-n)))) * sizeof(double));

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );
	lwork = ((int) work[0]);
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, A, lda, R, ldr );

	LAPACKE_dlaset_work( LAPACK_COL_MAJOR, 'A', m, n, 0.0e+00, 1.0e+00, Q, ldq );

//	dlamtsqr_( "L", "N", &m, &n, &n, &mb, &nb, A, &lda, T, &ldt, Q, &ldq, work, &lwork, &info );

	int k, kk, ctr, ii, i0, itmp;

	i0 = 0;

	k = n;

	if( ( mb <= k ) || ( ( mb >= m ) && ( mb >= n ) && ( mb >= k ) ) ) {

		dgemqrt_( "L", "N", &m, &n, &k, &nb, A, &lda, T, &ldt, Q, &ldq, work, &info);

	} else {

	kk = (m-k) % (mb-k);
	ctr = (m-k) / (mb-k);

//	printf(" kk = %d, ctr = %d\n", kk, ctr );

	if (kk > 0){

		ii=m-kk+1;

		dtpmqrt_( "L", "N", &kk, &n, &k, &i0, &nb, &(A[ii-1]), &lda, &(T[ ctr*k*ldt ]), &ldt, &(Q[0]), &ldq, &(Q[ii-1]), &ldq, work, &info );

	}
	else {

		ii=m+1;
	}

	for ( i = ii - (mb-k); i >= mb+1; i -= ( mb - k ) ){

//		printf(" i = %d, mb-k = %d, ctr = %d\n", i, mb-k, ctr );

		ctr = ctr - 1;

		itmp = mb - k;
		dtpmqrt_( "L", "N", &itmp, &n, &k, &i0, &nb, 
			&(A[i-1]), &lda, &(T[ctr*k*ldt]), &ldt, &(Q[0]), &ldq, &(Q[(i-1)]), &ldq, work, &info );
	}

	dgemqrt_( "L", "N", &mb, &n, &k, &nb, &(A[0]), &lda, &(T[0]), &ldt, &(Q[0]), &ldq, work, &info );

	}

	free( work );
	free( T );
*/



	int info;

 	T = (double *) malloc( ldt * n * ((int) ceil(( (double) (m-n) )/((double) (mb-n)))) * sizeof(double));

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );
	lwork = ((int) work[0]);
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, A, lda, R, ldr );

	LAPACKE_dlaset_work( LAPACK_COL_MAJOR, 'A', m, n, 0.0e+00, 1.0e+00, Q, ldq );

//	dlamtsqr_( "L", "N", &m, &n, &n, &mb, &nb, A, &lda, T, &ldt, Q, &ldq, work, &lwork, &info );

	int k, kk, ctr, ii, i0, itmp;
	int ib, kf;
	int jb, lb, iii;

	i0 = 0;

	k = n;

	if( ( mb <= k ) || ( ( mb >= m ) && ( mb >= n ) && ( mb >= k ) ) ) {

//		dgemqrt_( "L", "N", &m, &n, &k, &nb, A, &lda, T, &ldt, Q, &ldq, work, &info);

//		printf("===> here (1) \n");

		kf = ((k-1)/nb)*nb+1;

		for ( i = kf-1; i >= 0; i -= nb ){

			if ( nb < k - i ) ib = nb; else ib = k-i;

			LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', m-i, n, ib, &(A[ i + i*lda ]), lda, &(T[i*ldt]), ldt, &(Q[i]), ldq, work, n );

		}



	} else {

		kk = (m-k) % (mb-k);
		ctr = (m-k) / (mb-k);

//		printf(" kk = %d, ctr = %d\n", kk, ctr );

		if (kk > 0){

			ii=m-kk+1;

//			dtpmqrt_( "L", "N", &kk, &n, &k, &i0, &nb, &(A[ii-1]), &lda, &(T[ ctr*k*ldt ]), &ldt, &(Q[0]), &ldq, &(Q[ii-1]), &ldq, work, &info );
//			dtpmqrt_( "L", "N", &itmp, &n, &k, &i0, &nb, &(A[i-1]), &lda, &(T[ctr*k*ldt]), &ldt, &(Q[0]), &ldq, &(Q[(i-1)]), &ldq, work, &info );

//			SUBROUTINE DTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO )
//			L is 0

//			KF = ((K-1)/NB)*NB+1
			kf = ((k-1)/nb)*nb+1;

//			DO I = KF, 1, -NB
			for (iii  = kf; iii >= 1; iii -= nb ){

//				printf("===> here (11) \n");

//				IB = MIN( NB, K-I+1 )
				if ( nb < k - iii + 1 ) ib = nb; else ib = k-iii+1;

//				MB = MIN( M-L+I+IB-1, M )
//				if (  itmp+iii+ib-1 < itmp ) jb = itmp+iii+ib-1; else jb = itmp;
				if (  kk+iii+ib-1 < kk ) jb = kk+iii+ib-1; else jb = kk;

//				IF( I.GE.L ) THEN LB = 0 ELSE LB = MB-M+L-I+1 END IF
				if( iii >= 0 ) lb = 0; else lb = jb-kk-iii+1;

//				CALL DTPRFB( 'L', 'N', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB )
//				LAPACKE_dtprfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', jb, n, ib, lb, &(A[(i-1)+(iii-1)*lda]), lda, &(T[(iii-1)*ldt+ctr*k*ldt]), ldt, &(Q[(iii-1)]), ldq, &(Q[(i-1)]), ldq, work, ib );
				LAPACKE_dtprfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', jb, n, ib, lb, &(A[(ii-1)+(iii-1)*lda]), lda, &(T[(iii-1)*ldt+ctr*k*ldt]), ldt, &(Q[(iii-1)]), ldq, &(Q[(ii-1)]), ldq, work, ib );

//			END DO
			}






		}
		else {

			ii=m+1;
		}


		for ( i = ii - (mb-k); i >= mb+1; i -= ( mb - k ) ){

//			printf(" i = %d, mb-k = %d, ctr = %d\n", i, mb-k, ctr );

			ctr = ctr - 1;

			itmp = mb - k;
//			dtpmqrt_( "L", "N", &itmp, &n, &k, &i0, &nb, &(A[i-1]), &lda, &(T[ctr*k*ldt]), &ldt, &(Q[0]), &ldq, &(Q[(i-1)]), &ldq, work, &info );

//			SUBROUTINE DTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO )
//			L is 0

//			KF = ((K-1)/NB)*NB+1
			kf = ((k-1)/nb)*nb+1;

//			DO I = KF, 1, -NB
			for (iii  = kf; iii >= 1; iii -= nb ){

//				printf("===> here (11) \n");

//				IB = MIN( NB, K-I+1 )
				if ( nb < k - iii + 1 ) ib = nb; else ib = k-iii+1;

//				MB = MIN( M-L+I+IB-1, M )
				if (  itmp+iii+ib-1 < itmp ) jb = itmp+iii+ib-1; else jb = itmp;

//				IF( I.GE.L ) THEN LB = 0 ELSE LB = MB-M+L-I+1 END IF
				if( iii >= 0 ) lb = 0; else lb = jb-itmp-iii+1;

//				CALL DTPRFB( 'L', 'N', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB )
				LAPACKE_dtprfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', jb, n, ib, lb, &(A[(i-1)+(iii-1)*lda]), lda, &(T[(iii-1)*ldt+ctr*k*ldt]), ldt, &(Q[(iii-1)]), ldq, &(Q[(i-1)]), ldq, work, ib );

//			END DO
			}

		}

//		dgemqrt_( "L", "N", &mb, &n, &k, &nb, &(A[0]), &lda, &(T[0]), &ldt, &(Q[0]), &ldq, work, &info );

//		printf("===> here (2) \n");

		kf = ((k-1)/nb)*nb+1;

		for ( i = kf-1; i >= 0; i -= nb ){

			if ( nb < k - i ) ib = nb; else ib = k-i;

				LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', mb-i, n, ib, &(A[ i + i*lda ]), lda, &(T[i*ldt]), ldt, &(Q[i]), ldq, work, n );

		}
	}

	free( work );
	free( T );


//	perform = ((double) flops_lapack_org2r( m, n, k ) + (double) flops_lapack_geqr2( m, k ) ) / elapsed / 1.0e+9 ;
//	perform = 0.0e+00 ;


	printf("%6d %6d ", m, n);
//	printf("%16.8f %10.3f ", elapsed, perform);


	check_qq_orth( &orth, m, n, Q, ldq );
	printf(" %5.1e  ",orth); 

	check_qr_repres( &repres, m, n, Asave, lda, Q, ldq, R, ldr );
	printf(" %5.1e  ",repres); 

	printf("\n");		

	free( Asave );
	free( A );
	free( Q );
	free( R );

	return 0;

}

int check_qr_repres( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr ){

	double normA, *work;
	int ii, jj;

	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, A, lda, NULL );

	work  = (double *) malloc(m * n * sizeof(double));
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), R, ldr, work, m );
 	for(ii = 0; ii < m; ii++) for(jj = 0; jj < n; jj++) work[ ii+jj*m ] -= A[ ii+jj*lda ];
	(*norm_repres) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	free( work );

	(*norm_repres) = (*norm_repres) / normA;

	return 0;
}

int check_qq_orth( double *norm_orth_1, int m, int n, double *Q, int ldq  ){

	double *work;

	work  = (double *) malloc( n * n * sizeof(double));
	LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
	(*norm_orth_1) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );

	return 0;
}
