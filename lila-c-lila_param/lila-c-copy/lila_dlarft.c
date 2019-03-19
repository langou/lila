#include "lila.h"

int lila_dlarft( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau ){

	int vrtq;
	vrtq = lila_param[4];

	// vrtq = 0  ==> VRTQ
	// vrtq = 1  ==>  VRT
	// vrtq = 3  ==>    T
	if(( vrtq == 0 ) || ( vrtq == 1 ) || ( vrtq == 3 )){

		double *Tki, *Aii, *T0j, *Ajj;
		int vb, info, k, j, ml;

		k = i % mt;

		Tki = T + k + i*ldt;
		Aii = A + i + i*lda;

		ml = m - i;
		vb = mt - k; if ( vb > n ) vb = n;

//		double *Akk; 
//		double normv2; 
//		Akk=Aii;
//		for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); printf("  %f %f\n", tau[k], 2.0e+00 / normv2); Akk=Akk+1+lda; }
//		Akk=Aii;
//		for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Aii, lda, tau, Tki, ldt);

		j = i + vb;
		ml -= vb;

		Ajj = A + j + j*lda;
		T0j = T     + j*ldt;
		tau += vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

		while( vb != 0 ){

//		Akk=Ajj;
//		for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); printf("  %f %f\n", tau[k], 2.0e+00 / normv2); Akk=Akk+1+lda; }
//		Akk=Ajj;
//		for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

			info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Ajj, lda, tau, T0j, ldt);
		
			j += vb;
			ml -= vb;

			Ajj += ( vb + vb * lda ) ;
			T0j += ( vb*ldt );
			tau += vb;

			if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

		}

	}

	return 0;
}   
