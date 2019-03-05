#include "lila.h"

int lila_query_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S ){

	int lwork1;
	lwork1 = 0;

	if( lwork1 < n*(j-i)+n) lwork1 = n*(j-i)+n;

	return lwork1; 

}

