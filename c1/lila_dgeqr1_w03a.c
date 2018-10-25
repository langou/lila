#include "lila.h"

int lila_dgeqr1_w03a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Aii, *Qii, *Tii;
	int ml;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + (i % mt) + i*ldt;

	ml = m - i;

  	cblas_dcopy( ml, Aii, 1, Qii, 1);

  	LAPACKE_dlarfg( ml, Aii, Aii+1, 1, Tii );

	cblas_dscal( ml, ( 1.0e+00 / (*Aii) ), Qii, 1);

	return 0;
}
