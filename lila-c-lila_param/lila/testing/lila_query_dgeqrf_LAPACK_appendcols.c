#include "lila.h"

int lila_query_dgeqrf_LAPACK_appendcols( int m, int k, int n, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	int lwork1;

	lwork1 = 0;	
	lwork1 = n*k; if ( lwork < lwork1 ) lwork = lwork1;
	lwork1 = n*n; if ( lwork < lwork1 ) lwork = lwork1;

	return lwork;

}
