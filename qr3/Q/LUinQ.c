#include "lila.h"

int LUinQ( int n, double *A, int lda ){

	int info, n1, n2; 
	double *work, *U11, *U12, *U22, *A11, *A12, *A21, *A22;

	if ( n <= 1 ){



	} else {
 
		// This is the cheat
		work = (double *) malloc( n * n * sizeof(double));
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, A, lda, work, n ); 
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), A, lda);
		

		n1 = n/3; n2 = n-n1; // doing n/3 for checking to make blocks not equal. If not changed through recursion you will get a segfault

		A11 = A;
		U11 = work;

		U12 = work + n1*n;
		A12 = A + n1*lda;
		A21 = A + n1; 

		U22 = work + n1 + n1*n;  
		A22 = A + n1 + n1*lda; 




//////////   Not sure which trmm (potentially gemms) can be removed from recursion. Also, we need 
//////////   the original T*V1^T to keep being put in the top of A through this routine. 

		// The computation leaves the square part of Q dense. So we break into 4
		// Q22
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n2, n2, (-1.0e+00), U22, n, A22, lda );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n2, n2, n1, (-1.0e+00), A21, lda, U12, n, (+1.0e+00), A22, lda );
		// Q12
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n1, n2, n1, (-1.0e+00), A11, lda, U12, n, (+0.0e+00), A12, lda );
		// Q11
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n1, (-1.0e+00), U11, n, A11, lda );
		// Q21
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n2, n1, (-1.0e+00), U11, n, A21, lda );


	//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, work, n, A, lda ); 
	//	info = LUinQ( n1, A11, lda );
	//	info = LUinQ( n2, A22, lda );




		free( work );
	}



	return 0;

}


