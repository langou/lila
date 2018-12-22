#include "lila.h"

int dgmres( int n, double *b, double *x, int m, int max_it, double tol ){

	double *H, *cs, *sn, *s;
	double *Q;

	double *Mlb, *Mlr, *Mrx, *AMrx, *MlAMrx;

	double *arnoldi_res, *true_res;

	double *H_;

	double bnrm2;
	int i, info, iter;

	int mgs_lvl1;
	int check;

	mgs_lvl1 = 1;
	check = 1;

	H  = (double *) malloc( (m+1) * m * sizeof(double) );
	cs = (double *) malloc( m * sizeof(double) );
	sn = (double *) malloc( m * sizeof(double) );
	s  = (double *) malloc( (m+1) * sizeof(double) );

	arnoldi_res = (double *) malloc( (m+1) * sizeof(double) );
	true_res = (double *) malloc( (m+1) * sizeof(double) );

	Mlr    = (double *) malloc( n * sizeof(double) );
	Mrx    = (double *) malloc( n * sizeof(double) );
	AMrx   = (double *) malloc( n * sizeof(double) );
	MlAMrx = (double *) malloc( n * sizeof(double) );
	Mlb    = (double *) malloc( n * sizeof(double) );

if (mgs_lvl1){
	Q = (double *) malloc( n * (m+1) * sizeof(double) );
	if ( check ) H_ = (double *) malloc( (m+1) * m * sizeof(double) );
}

	iter = 1;
//
//	bnrm2 = norm(Ml*b,2);
	info = matvec_Ml( n, Mlb, b );
	bnrm2 = cblas_dnrm2( n, Mlb, 1);
//
//	r = Ml * ( b-A*Mr*x);
	info = matvec_Mr( n, Mrx, x );
	info = matvec_A ( n, AMrx, Mrx );
	info = matvec_A ( n, MlAMrx, AMrx );
 	for(i = 0; i < n; i++) Mlr[i] = Mlb[i] - MlAMrx[i];

if (mgs_lvl1){
 	for(i = 0; i < n; i++) Q[i] = Mlr[i];
	s[0] = cblas_dnrm2( n, Q, 1);
 	for(i = 0; i < n; i++) Q[i] /= s[0];
}

	arnoldi_res[0] = fabs( s[0] ) / bnrm2;
	true_res[0] = fabs( s[0] ) / bnrm2;


	printf("\n");
	printf("%e %e", arnoldi_res[0], true_res[0] );
	printf("\n");

	

	




if (mgs_lvl1){
	free( Q );
	if ( check ) free( H_ );
}

	free( arnoldi_res );
	free( true_res );

	free( Mlb );
	free( MlAMrx );
	free( AMrx );
	free( Mrx );
	free( Mlr );

	free(s);
	free(sn);
	free(cs);
	free(H);

	return 0;

}
