#include "lila.h"

int lila_dgeqrf_qr2( int *lila_param,  int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info;

	if( lila_param[4] == 0 ){

		if     ( lila_param[2] == 0 ){	
			info = lila_dgeqr2_w03_l ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}
		else if( lila_param[2] == 1 ){
			info = lila_dgeqr2_w03_3 ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}
		else if( lila_param[2] == 2 ){
			info = lila_dgeqr2_w03_hr( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}

	} else if( lila_param[4] == 1 ){

		if     ( lila_param[2] == 0 ){	
			info = lila_dgeqr2_v03_l ( m, n, i, mt, A, lda, T, ldt, work, lwork );
		}
		else if( lila_param[2] == 1 ){
			info = lila_dgeqr2_v03_3 ( m, n, i, mt, A, lda, T, ldt, work, lwork );
		}
		else if( lila_param[2] == 2 ){
		//	info = lila_dgeqr2_v03_hr( m, n, i, mt, A, lda, T, ldt, work, lwork );
		}

	} else if( lila_param[4] == 2 ){

		if     ( lila_param[2] == 0 ){	
			info = lila_dgeqr2_q03_l ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}
		else if( lila_param[2] == 1 ){
			info = lila_dgeqr2_q03_3 ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}
		else if( lila_param[2] == 2 ){
		//	info = lila_dgeqr2_q03_hr( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}

	}

	return 0;
}
