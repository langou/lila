#include "lila.h"

int lila_wsq_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S ){
	int lwork1;
	lwork1 = 0;

	//if( lwork1 < (n+i)*(j-i)) lwork1 = (n+i)*(j-i);
	//if( lwork1 < n*(j-i)) lwork1 = n*(j-i);
	//if( lwork1 < (n+i)*n) lwork1 = (n+i)*n;

	if( lwork1 < n*(j-i)+n) lwork1 = n*(j-i)+n;// I think this is the right size we need for work
						   // but testing, it doesn't seem to change anything by commenting them all out.
						   // a bigger workspace must be needed in dormqrf and dormqrbz
	return lwork1; 

}

