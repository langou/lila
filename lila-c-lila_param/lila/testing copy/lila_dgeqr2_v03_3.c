#include "lila.h"

int lila_dgeqr2_v03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *work1, *Aii, *Tki;
	int k, ml, info; 

	work1 = work;	

	ml = m - i;
	k  =  i%mt;

	Aii = A + i*lda + i;
	Tki = T + i*ldt + k;

  	info = dgeqr3( ml, n, Aii, lda, Tki, ldt );

	return 0;

}
