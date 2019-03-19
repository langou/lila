#include "lila.h"

int lila_dgeqr2( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, panel;
	panel = lila_param[2];

	if     ( panel == 0 ){	
		info = lila_dgeqr2_l ( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}
	else if( panel == 1 ){
		info = lila_dgeqr2_3 ( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}
	else if( panel == 2 ){
		info = lila_dgeqr2_hr( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	return 0;
}
