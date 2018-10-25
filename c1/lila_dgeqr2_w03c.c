#include "lila.h"

int lila_dgeqr2_w03c( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii, *TTTii;
	int ml, ii, accum;


	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	TTTii = TTT + (i % mt) + i*llldddttt;

	ml = m - i;

  	cblas_dcopy( ml, Aii, 1, Qii, 1);

  	LAPACKE_dlarfg( ml, Aii, Aii+1, 1, TTTii );

	cblas_dscal( ml, ( 1.0e+00 / (*Aii) ), Qii, 1);

//	as long as we want to maintain T
	Tii = T + i*ldt + i;
	(*Tii) = (*TTTii);

	return 0;
}
