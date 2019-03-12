#include "lila.h"

int lila_dV2tau_w03( int m, int n, int i, int mt, double *A, int lda, double *tau ){

	double *Akk, *Ajj, normv2; 
	int ml, vb, k, j;

	// This script works if mt == n. I have not tested it on another mt

	ml = m-i;
	vb  = mt - (i%mt); if ( vb > n ) vb = n;
	Akk = A + i + i*lda;

	for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

	j    = i + vb;
	tau += vb;
	Ajj = A + (i+vb) + (i+vb)*lda;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		Akk = Ajj;
		for( k = 0; k < vb; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

		j   += vb;
		tau += vb;
		Ajj += vb*(1+lda);

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}   
