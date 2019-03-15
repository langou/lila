#include "lila.h"

int lila_dgeqr2_v03( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	int info, panel;
	panel = lila_param[2];

	if( panel == 0 ){	
		info = lila_dgeqr2_v03_l ( m, n, i, mt, A, lda, T, ldt, work, lwork );
	}
	else if( panel == 1 ){
		info = lila_dgeqr2_v03_3 ( m, n, i, mt, A, lda, T, ldt, work, lwork );
	}
	else if( panel == 2){
		//info = lila_dgeqr2_w03_hr( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	return 0;
}
