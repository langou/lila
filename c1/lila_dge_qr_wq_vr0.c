#include <stdio.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int lila_dge_qr_wq_vL0( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

int lila_dge_qr_wq_vr0( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq ){

	int j, info, lwork, lwork_max; 
	double *work=NULL;
	double *tau=NULL;

	tau = (double *) malloc( n * sizeof(double));

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork );
	printf("lwork GEQRF = %d\n",((int) work[0]));
	lwork_max = ((int) work[0]);
	lwork = -1;
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, A, lda, tau, work, lwork );
	printf("lwork ORGQR = %d\n",((int) work[0]));
	if( ((int) work[0]) > lwork_max ) lwork_max = ((int) work[0]); 
	lwork = lwork_max;
	printf("lwork       = %d\n",lwork);

	work = (double *) malloc( lwork * sizeof(double));

	int n1, n2;

	n1 = n / 2;
	n2 = n - n1;

  	info = lila_dge_qr_wq_vL0( m, n1, 0, A, lda, T, ldt, Q, ldq, work, lwork );

  	info = lila_dge_qr_ormqrf_vL0( m, n2, n1, 0, n1, mt, A, lda, T, ldt );

	for(j = 0; j < n1; j++) tau[j] = T[j+j*ldt];
	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', m, n2, n1, A, lda, tau, &(A[n1*lda]), lda, work, lwork );

  	info = lila_dge_qr_wq_vL0( m, n2, n1, A, lda, T, ldt, Q, ldq, work, lwork );

//	for(j = n1; j < n; j++) tau[j] = T[j+j*ldt];
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );
//	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n1, n2, (0e+00), (0e+00), &(Q[n1*ldq]), ldq );
	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'N', m, n2, n1, A, lda, tau, &(Q[n1*ldq]), ldq, work, lwork );

	free( work );
	free( tau );

	return 0;

}
