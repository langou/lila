#include "lila.h"

int lila_dgeqr1_w03b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	double *Aii, *Tii, *Qii, *TTTii;
	double norma, normv;
	double *Aki, *Qki;
	int kk;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	TTTii = TTT + (i % mt) + i*llldddttt;

	for ( kk = i, Aki = Aii, Qki = Qii; kk < m; kk++, Aki++, Qki++ ) (*Qki) = (*Aki);

	normv = 0e+00; for ( kk = i+1, Aki = Aii+1; kk < m; kk++, Aki++ ) normv += (*Aki)*(*Aki);
	norma = normv; norma += (*Aii)*(*Aii); norma = sqrt( norma );

	if( (*Aii) > 0 ) (*Aii) += norma; else  (*Aii) -= norma;

	(*TTTii) = 2.0e+00 / ( 1.0e+00 + ( ( normv / (*Aii) ) / (*Aii) ) );

	for ( kk = i+1, Aki = Aii+1; kk < m; kk++, Aki++ ) (*Aki) /= (*Aii);

	if( (*Aii) > 0 ) (*Aii) = -norma; else  (*Aii) = +norma;

	for ( kk = i, Qki = Qii; kk < m; kk++, Qki++ ) (*Qki) /= (*Aii);

//	as long as we need to maintain T as well as TTT, these two lines will go away at some point
	Tii = T + i*ldt + i;
	(*Tii) = (*TTTii);

	return 0;
}
