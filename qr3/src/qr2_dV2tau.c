#include "qr2.h"

int qr2_dV2tau( int m, int n, double *A, int lda, double *tau ){

	double *Akk, normv2; 
	int k;

	Akk = A;

	for( k = 0; k < n; k++){ normv2=1+cblas_ddot(m-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

	return 0;
}   
