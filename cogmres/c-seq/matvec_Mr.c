#include "lila.h"

int matvec_Mr( int n, double *y, double *x ){

	int i;

 	for(i = 0; i < n; i++){
		y[i] = 1.e-1 * ( (double) i+1 ) * x[i];
	}

	return 0;
}