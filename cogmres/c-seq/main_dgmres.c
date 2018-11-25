#include "lila.h"

extern int matvec_A( int n, double *y, double *x );
extern int matvec_Ml( int n, double *y, double *x );

int main(int argc, char ** argv) {

	int i, n, info;
	double *x, *y, *b;
	double elapsed_refL;
	struct timeval tp;

	n = 20;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("n = %4d; ",n);

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


///////

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
//	printf("time = %f;", elapsed_refL );
//	printf("perf = %f; \n", 2.0e+00*((double) m)*((double) n)*((double) k)/elapsed_refL*1e-9 );

	free(x);
	free(b);

	return 0;
}
