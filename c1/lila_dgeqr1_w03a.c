#include "lila.h"

int lila_dgeqr1_w03a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	double *Aii, *Tii, *Qii, *TTTii;
	int ml;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	TTTii = TTT + (i % mt) + i*llldddttt;

	ml = m - i;

  	cblas_dcopy( ml, Aii, 1, Qii, 1);

  	LAPACKE_dlarfg( ml, Aii, Aii+1, 1, TTTii );

	cblas_dscal( ml, ( 1.0e+00 / (*Aii) ), Qii, 1);

//	as long as we need to maintain T as well as TTT, these two lines will go away at some point
	Tii = T + i*ldt + i;
	(*Tii) = (*TTTii);

	return 0;
}
