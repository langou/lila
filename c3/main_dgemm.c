#include "lila.h"

int main(int argc, char ** argv) {

	int i; double *tmp;
	int lda, ldb, ldc, m, n, k;
	double *A, *B, *C;
	double elapsed_refL;
	struct timeval tp;

	srand(0);

	m = 5;
	n = 10;
	k = 20;
	lda = 20;
	ldb = 30;
	ldc = 30;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-k") == 0) {
			k  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldb") == 0) {
			ldb  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldc") == 0) {
			ldc  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if (lda < m) lda=m;
	if (ldb < k) ldb=k;
	if (ldc < m) ldc=m;

	printf("m = %4d; ",m);
	printf("n = %4d; ",n);
	printf("k = %4d; ",n);
	printf("lda = %4d; ",lda);
	printf("ldb = %4d; ",ldb);
	printf("ldc = %4d; ",ldc);

	B  = (double *) malloc(ldb * n * sizeof(double));
	C = (double *) malloc(ldc * n * sizeof(double));
	A  = (double *) malloc(lda * k * sizeof(double));

 	for(i = 0, tmp=B; i < ldb * n; i++, tmp++){
		*tmp = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	}
 	for(i = 0, tmp=A; i < lda * k; i++, tmp++){
		*tmp = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	}
 	for(i = 0, tmp=C; i < ldc * n; i++, tmp++){
		*tmp = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	}

/*	call reference LAPACK trsm */

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
		m, n, k, 1.0e+00, A, lda, B, ldb, 1.0e+00, C, ldc);
	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	printf("time = %f;", elapsed_refL );
	printf("perf = %f; \n", 2.0e+00*((double) m)*((double) n)*((double) k)/elapsed_refL*1e-9 );

	free(A);
	free(C);
	free(B);

	return 0;
}
