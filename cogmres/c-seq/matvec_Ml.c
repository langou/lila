#include "lila.h"

int matvec_Ml( int n, double *y, double *x ){

	int i;

 	for(i = 0; i < n; i++){
		y[i] = (1.0e+00 + 1.e-1 * ( (double) i+1 )) * x[i];
	}

	return 0;
}
