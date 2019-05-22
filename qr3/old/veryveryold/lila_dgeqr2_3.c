#include "lila.h"

int lila_dgeqr2_3( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vrtq;
	vrtq = lila_param[4];

	// vrtq == 0  ==> VRTQ
	if( vrtq == 0 ){

		double *work1, *Aii, *Qii, *Tki;
		int j, k, ml, info; 

		work1 = work;	

		ml = m - i;
		k  =  i%mt;

		Aii = A + i*lda + i;
		Qii = Q + i*ldq + i;
		Tki = T + i*ldt + k;

	  	info = dgeqr3( ml, n, Aii, lda, Tki, ldt );

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Qii, ldq);
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work1, n);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Tki, ldt, work1, n );

//		note that this is a triangle times triangle so we could have dtrtrmm() and dived # of FLOPS by 3x
		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Aii, lda, work1, n );
	
//		note that the top n-by-n part is a lower times a upper and can be done dlauum
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (-1.0e+00), work1, n, Qii, ldq );

	 	for(j = 0; j < n; j++) Qii[ j + ldq * j ] = 1.00e+00 + Qii[ j + ldq * j ];

	}

	// vrtq == 1  ==> VRT
	if( vrtq == 1 ){

		double *work1, *Aii, *Tki;
		int k, ml, info; 

		work1 = work;	

		ml = m - i;
		k  =  i%mt;

		Aii = A + i*lda + i;
		Tki = T + i*ldt + k;

	  	info = dgeqr3( ml, n, Aii, lda, Tki, ldt );

	}

	// vrtq == 2  ==> Q
	if( vrtq == 2 ){

		double *work1, *Aii, *Qii, *Tki;
		int j, k, ml, info; 

		work1 = work;	

		ml = m - i;
		k  =  i%mt;

		Aii = A + i*lda + i;
		Qii = Q + i*ldq + i;
		Tki = T + i*ldt + k;

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Qii, ldq);
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work1, n);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Tki, ldt, work1, n );

//		note that this is a triangle times triangle so we could have dtrtrmm() and dived # of FLOPS by 3x
		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Aii, lda, work1, n );

//		note that the top n-by-n part is a lower times a upper and can be done dlauum
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (-1.0e+00), work1, n, Qii, ldq );

	 	for(j = 0; j < n; j++) Qii[ j + ldq * j ] = 1.00e+00 + Qii[ j + ldq * j ];

	}

	// THIS ONE NEEDS FIXING 
	// vrtq == 3  ==> T
	if( vrtq == 3 ){

		double *Aii, *Akk, *Tki, *tau=NULL, normv2;
		int j, k, ml, info, kk; 

		ml = m - i;

		tau = work;
		Aii = A + i*lda + i;
		Tki = T + i*ldt + (i%mt);
		Akk = Aii;
		for( kk = 0; kk < n; kk++){ normv2=1+cblas_ddot(ml-kk-1,Akk+1,1,Akk+1,1); tau[kk] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

		for( k = 0; k < n; k++){

	  		info = dgeqr3( ml, 1, tau, 1, Tki, ldt );
			ml -= 1;
			j  += 1;
			tau += 1;
			Tki = T + j*ldt + (i%mt);

		}

	}

	return 0;

}
