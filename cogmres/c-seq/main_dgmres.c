#include "lila.h"

int main(int argc, char ** argv) {

	int i, n, info, m, max_it;
	double *x, *y, *b;
	double elapsed_refL;
	double tol;
	struct timeval tp;

	n = 20;
	m = 5;
	max_it = 10;
	tol = 1.0e-6;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-tol") == 0) {
			tol  = atof( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-max_it") == 0) {
			max_it  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("n = %4d; tol = %5.0e; m = %2d; max_it = %3d; ", n, tol, m, max_it );

	b  = (double *) malloc(n * sizeof(double));
	x  = (double *) malloc(n * sizeof(double));

 	for(i = 0; i < n; i++){
		x[i] = 1.00e+00;
	}
 	for(i = 0; i < n; i++){
		b[i] = ((double) (n))/(2.0e+00) - ((double) (i+1));
	}

	y  = (double *) malloc(n * sizeof(double));
//	info = matvec_A( n, y, b );
//	info = matvec_Ml( n, y, b );
	info = matvec_Mr( n, y, b );
	printf("y=[ "); for(i = 0; i < n-1; i++) printf("%f, ", y[i]); printf("%f ];\n",y[n-1]); 
	free(y);

/*	call reference LAPACK trsm */

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);


	info = dgmres( n, b, x, m, max_it, tol );

///////

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
//	printf("time = %f;", elapsed_refL );
//	printf("perf = %f; \n", 2.0e+00*((double) m)*((double) n)*((double) k)/elapsed_refL*1e-9 );

	free(x);
	free(b);

	return 0;
}
