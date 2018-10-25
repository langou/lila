#include "lila.h"

int lila_dge_qr_ormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Qij, *Tii;
	int vb, jj, not_done, ldwork, itlo, ml;

	vb = ((i+k) % mt); if (vb == 0) vb = mt; if ( vb > k ) vb = k;

//	This is the logic that follows the code in matlab
	Aii = A + (i+k-vb) + (i+k-vb)*lda; 
	Qij = Q + (i+k-vb) + j*ldq;
	Tii = T + (i % mt) + (i+k-vb)*ldt;

	itlo = ((i+k) % mt);

	ldwork = k;

	ml = m - i - k + vb;
	printf("i = %d, k = %d, vb = %d, itlo = %d, ml = %d,\n",i,k,vb,itlo,ml);

	not_done = 1;

	jj = 1;

//	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', k, n, (0e+00), (0e+00), Qij, ldq );

//	The first go works, tracking within the while is where I am currently struggling	
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, m-i-vb, (1.0e+00), Aii+vb, lda, Qij+vb, lda, (0.0e+00), work, ldwork );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', k, vb, work, ldwork, Qij, ldq ); // k, vb    or    vb, vb   ?
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (-1.0e+00), Aii, lda, Qij, ldq );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-i-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );

	while ( not_done == 1 ){

		if ( jj + vb - 1 == k ){ 

			not_done = 0; 

		} else if ( jj + vb - 1 > k ){

			printf("defensive programming, we should never have been there, abort\n");

		} else{
			printf("\nTHIS IS IN THE WHILE\n");	

			jj += vb;

			if ( ( jj + mt - 1 ) <= k ) vb = mt; else vb = k - jj + 1; 

//			itlo = ( (itlo-vb) % mt ); if (itlo == 0) itlo = mt;
			itlo = ( abs((itlo-vb)) % mt ); if (itlo == 0) itlo = mt; // does Matlab automatically make mod() positive?
			ml += vb;

			printf("itlo = %d, vb = %d, k = %d\n",itlo,vb,k);

			Aii = Aii - vb - vb*lda ;
			Qij = Qij - vb;
			Tii = Tii - itlo - vb*ldt;
//			Tii = Tii - ((i-vb)%mt) - vb*ldt;

			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml, (1.0e+00), Aii+vb, lda, Qij+vb, lda, (0.0e+00), work, ldwork );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
			LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', k, vb, work, ldwork, Qij, ldq ); // k, vb    or    vb, vb   ?
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (-1.0e+00), Aii, lda, Qij, ldq );
			cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );
		}

	}

	return 0;


/*
	vb = ((i+k) % mt); if (vb == 0) vb = mt; if ( vb > k ) vb = k;

//	This is the logic that follows the code in matlab
	Aii = A + (i+k-vb) + (i+k-vb)*lda; 
	Qij = Q + (i+k-vb) + j*ldq;
	Tii = T + (i % mt) + (i+k-vb)*ldt;

	itlo = (i % mt);

	ldwork = k;

	ml = m - vb;

	not_done = 1;

	jj = 1;

	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml, (1.0e+00), Aii+vb, lda, Qij+vb, lda, (0.0e+00), work, ldwork );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', k, vb, work, ldwork, Qij, ldq ); 
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (-1.0e+00), Aii, lda, Qij, ldq );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );

	while ( not_done == 1 ){

		if ( jj + vb - 1 == k ){ 

			not_done = 0; 

		} else if ( jj + vb - 1 > k ){

			printf("defensive programming, we should never have been there, abort\n");

		} else{
			printf("\nTHIS IS IN THE WHILE\n");	

			jj += vb;

			if ( ( jj + mt - 1 ) <= k ) vb = mt; else vb = k - jj + 1; 

//			itlo = ( (itlo-vb) % mt ); if (itlo == 0) itlo = mt;
			itlo = ( abs((itlo-vb)) % mt ); if (itlo == 0) itlo = mt; // does Matlab automatically make mod() positive?
			ml -= vb;

			Aii = Aii + vb + vb*lda ;
			Qij = Qij + vb;
			Tii = Tii + itlo + vb*ldt;

			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml, (1.0e+00), Aii+vb, lda, Qij+vb, lda, (0.0e+00), work, ldwork );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
			LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', k, vb, work, ldwork, Qij, ldq );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (-1.0e+00), Aii, lda, Qij, ldq );
			cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );
		}

	}
	return 0;
*/

}
