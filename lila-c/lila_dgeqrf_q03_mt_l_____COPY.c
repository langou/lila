#include "lila.h"

int lila_dgeqrf_q03_mt_l( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

printf("\n");

	double *tau=NULL, *Aii, *Qii, *Tki;
	int ml, info, k, vb, j, l, lwork1; 

	tau = work;
	work = work + n;
	lwork1 = lwork-n;


	//info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau);

	k = i % mt;

	Tki = T + k + i*ldt;
	Aii = A + i + i*lda;
	Qii = Q + i*ldq + i;

	ml = m - i;

	vb = mt - k; if ( vb > n ) vb = n;

	double *Akk; 
	double *Ajj; 
	double *T0j; 
	double *Tll; 
	double normv2; 

	Akk=Aii;
	Tll=Tki;
	for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); printf("  %f %f\n", *Tll, 2.0e+00 / normv2); Akk=Akk+1+lda; Tll=Tll+1+ldt; }
	Akk=Aii;
	Tll=Tki;
	for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; Tll=Tll+1+ldt; }

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Aii, lda, tau, Tki, ldt);

	j = i + vb;
	ml -= vb;

	Ajj = A + j + j*lda;
	T0j = T     + j*ldt;
	tau += vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){


	Akk=Ajj;
	Tll=T0j;
	for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); printf("  %f %f\n", *Tll, 2.0e+00 / normv2); Akk=Akk+1+lda; Tll=Tll+1+ldt; }
	Akk=Ajj;
	Tll=T0j;
	for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; Tll=Tll+1+ldt; }



		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Ajj, lda, tau, T0j, ldt);
		
		j += vb;
		ml -= vb;

		Ajj += ( vb + vb * lda ) ;
		T0j += ( vb*ldt );
		tau += vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}


	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

	return 0;

}
