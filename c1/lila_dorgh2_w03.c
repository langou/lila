#include "lila.h"

int lila_dorgh2_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){

	double *Tjj, *Ajj, *T0j;
	double *Tkk, *Q0j;
	int *Sj;
	int k, i1, ml, vb;

	vb = mt - (j % mt ); if ( vb > n ) vb = n;
	printf(" n = %d, (j%%mt) = %d, i = %d, j = %d, vb = %d \n",n, j%mt, i, j, vb);

	ml = m-i;

//	Tjj = T + j + j*ldt;
	Tjj = T + (j % mt) + j*ldt;
	T0j = T + j*ldt;

	Ajj = A + j + j*lda;
	Q0j = Q + j*ldq;
	Sj = S + j;

	int not_done, jj;
	not_done = 1;

	jj = 1;

	while( not_done == 1 ){

		for( k = 0; k < vb; k++ ){
			Tkk = Tjj + k + k*ldt;

			if ( fabs( 1.0e+00 - (*Tkk) ) < fabs( 1.0e+00 + (*Tkk) ) ){
				for( i1 = 0; i1 < vb; i1++) Tjj[ i1 + k*ldt ] = - Tjj[ i1 + k*ldt];
	 	 		Sj[ k ] = - 1;		
			} else {
				Sj[ k ] = 1;		
			}
			(*Tkk) = (*Tkk) - 1.0e+00;
			for( i1 = k+1; i1 < vb; i1++ ) Tjj[ i1 + k*ldt ]  = Tjj[ i1 + k*ldt ] / (*Tkk);
			cblas_dger( CblasColMajor, vb-k-1, vb-k-1, (-1.0e+00), Tkk+1, 1, Tkk+ldt, ldt, Tkk+1+ldt, ldt );

		}

		for( k = 0; k < vb; k++ ){
			if ( Sj[ k ] == -1 ){
				for( i1 = vb; i1 < m-j; i1++) Ajj[ i1 + k*lda ] = - Ajj[ i1 + k*lda ];
			}
		}

		cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml-j-vb, vb, 1.0e+00, Tjj, ldt, Ajj + vb, lda );

		for( k = 0; k < vb; k++ ){
			if ( Sj[ k ] == -1 ){
				for( i1 = 0; i1 < j-1; i1++) T0j[ i1 + k*ldt ] = - T0j[ i1 + k*ldt ];
			}
		}
		for( k = 0; k < vb; k++ ){
			if ( Sj[ k ] == -1 ){
				for( i1 = k; i1 < vb; i1++) Ajj[ k + i1*lda ] = - Ajj[ k + i1*lda ];
			}
		}
		for( k = 0; k < vb; k++ ){
			if ( Sj[ k ] == -1 ){
				for( i1 = 0; i1 < m-i; i1++) Q0j[ i1 + k*ldq ] = - Q0j[ i1 + k*ldq ];
			}
		}
		for( k = 0; k < vb; k++ ){
			for( i1 = 0; i1 < k; i1++) Ajj[ k + i1*lda ] = Tjj[ k + i1*ldt ];
		}
		for( k = 0; k < vb; k++ ){
			for( i1 = 0; i1 < k; i1++) Tjj[ k + i1*ldt ] = 0.0e+00;
		}

		cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, (-1.0e+00), Ajj, lda, Tjj, ldt );




		if ( jj + vb - 1 == n ) {

			not_done = 0; 

		} else if ( jj + vb - 1 > n ) {

			printf("defensive programming, we should never have been there, abort\n"); return 0;

		} else {

			if(jj == 1) Tjj = Tjj - (j%mt) + vb*ldt; 

			ml -= vb;
			jj += vb;

			Ajj += vb * ( lda + 1 );
			Sj += vb;

			if(jj != 1+vb) Tjj = Tjj + vb*ldt; 

			if ( ( jj + mt - 1 ) <= n ) vb = mt; else vb = n - jj + 1;

		}

	}

	return 0;

}
