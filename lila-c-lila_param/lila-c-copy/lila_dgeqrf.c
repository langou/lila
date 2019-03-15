int lila_dgeqrf( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	char mode;
	mode = lila_param[0];
	
	if ( mode == 'l' ){
		int n_lvl, *nb_lvl, j;
		n_lvl  = lila_param[ 3 ];
		nb_lvl = (int *) malloc(n_lvl * sizeof(int));
		for(j = 0; j < n_lvl; j++) nb_lvl[j] = lila_param[5+j];
		lila_dgeqrf_vrtq_levelx( lila_param, n_lvl, 0, nb_lvl, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );		
	}
	if ( mode == 'r' ){
		lila_dgeqrf_vrtq_recursive( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );		
	}

	return 0;
}

