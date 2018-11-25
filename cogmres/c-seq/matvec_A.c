#include "lila.h"

int matvec_A( int n, double *y, double *x ){

	int i;

	y[0] = 1.00e+01 * x[0] + 1.00e+00 * x[1];
 	for(i = 1; i < n-1; i++){
		y[i] = 2.00e+00 * x[i-1] + 1.00e+01 * x[i] + 1.00e+00 * x[i+1];
	}
	y[n-1] = 2.00e+00 * x[n-2] + 1.00e+01 * x[n-1];

	return 0;
}

// perform y <- A*x where A is 
//
// A = [ 10  1  0  0  ..   0  0
//        2 10  1  0  ..   0  0
//        0  2 10  1  ..   0  0
//        :  :  :  :   .   :  :
//        0  0  0  0  ..  10  1 
//        0  0  0  0  ..   2 10 ];
