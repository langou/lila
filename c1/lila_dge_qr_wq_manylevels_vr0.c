#include "lila.h"

int lila_dge_qr_wq_manylevels_vr0( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq ){

	int vb, info, lwork, lwork_max; 
	double *work=NULL;
	int nb;

	nb = nb_lvl[ i_lvl ];

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, NULL, work, lwork );
	printf("lwork GEQRF = %d\n",((int) work[0]));
	lwork_max = ((int) work[0]);
	lwork = -1;
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, A, lda, NULL, work, lwork );
	printf("lwork ORGQR = %d\n",((int) work[0]));
	if( ((int) work[0]) > lwork_max ) lwork_max = ((int) work[0]); 
	lwork = lwork_max;
	printf("lwork       = %d\n",lwork);

	work = (double *) malloc( lwork * sizeof(double));

	if ( nb > n ) vb = n; else vb = nb; 

	printf("vb (level 0) = %d\n", vb );

	info = lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0( n_lvl, 1, nb_lvl, m, vb, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//	info = lila_dge_qr_wq_vL0( m, vb, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	i = vb;
	
	while ( i < n ) {

	if ( i+nb > n ) vb = n-i; else vb = nb; 

	printf("vb (level 0) = %d\n", vb );

	info = lila_dge_qr_ormqrf_vL0( m, vb, i, 0, i, mt, A, lda, T, ldt, work, lwork );

	info = lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0( n_lvl, 1, nb_lvl, m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//	info = lila_dge_qr_wq_vL0( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	//info = lila_dge_qr_larft_connect_vL0( m, n, k, i, j, mt, A, lda, T, ldt, work, lwork );
	info = lila_dge_qr_ormqrbz_vL0( m, vb, i, 0, i, mt, A, lda, Q, ldq, T, ldt, work, lwork );
 
	i += vb;
	}

	free( work );

	return 0;

}
