#include "lila.h"

int lila_dgeqr2_w03c( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii, *TTTii;
	int ml, ii, accum;

	tau = (double *) malloc( n * sizeof(double));

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + i*ldt + i;
	TTTii = TTT + (i % mt) + i*llldddttt;

	ml = m - i;

	double normA;
	double ttt;

  	normA = cblas_dnrm2( ml, Aii, 1);

  	cblas_dcopy( ml, Aii, 1, Qii, 1);

//  	cblas_dscal( ml, ( 1.0e+00 / normA ), Qii, 1);

  	LAPACKE_dlarfg( ml, Qii, Qii+1, 1, &ttt );

	info = lila_dgeqr2_w03b( m, 1, i, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork ); // This function doesn't have TTT llldddttt in the interface

printf("===> %f ", ttt);
printf("|| %f || %f\n", *Tii, *TTTii);

	return 0;


//// LARFG  /////////////////

	// compute the norm
	for (ii = i; ii < i+n; ii++){
		accum += Aii[ii*lda]*Aii[ii*lda];	
	}

	if( Aii[i] > 0 ) Aii[i] += accum; else Aii[i] -= accum;

	for (ii = 1; ii < ml; ii++){
		Aii[ii*lda] = Aii[ii*lda] / Aii[i];
	}

	if( Aii[i] > 0 ) Aii[i] = -accum; else Aii[i] = accum;

//// LARFL  /////////////////


	
	free( tau );

	return 0;

}


/*

MATLAB CODE FOR LAPACK_GEQR2
	the loop of j goes from i to i+n

	[ A(j:m,j) ] = lapack_larfg( A(j:m,j) );
		m = size(a,1);
		norma = norm( a, 2);
		if ( a(1) > 0 ) a(1) = a(1) + norma; else a(1) = a(1) - norma; end
		a(2:m,1) = a(2:m) / a(1);
		if ( a(1) > 0 ) a(1) = - norma; else a(1) = norma; end
      
	[ A(j:m,j+1:ihi) ] = lapack_larfL( A(j:m,j), A(j:m,j+1:ihi) );
		[m,n]=size(B);
		tau = 2 / ( 1 + norm(a(2:m,1))^2 );
		tmp(1,1:n) = B(1,1:n) + a(2:m,1)' * B(2:m,1:n);
		tmp(1,1:n) = tau * tmp(1,1:n);
		B(1,1:n) = B(1,1:n) - tmp(1,1:n);
		B(2:m,1:n) = B(2:m,1:n) - a(2:m,1) * tmp(1,1:n);


MATLAB CODE FOR LILA_ORGQRF_V05_W00
	loop of j = i+n:-1:i

	Q(i:m,i:i+n-1) = eye( size ( Q(i:m,i:i+n-1) ) );
	for j = i+n-1:-1:i,
	      V = tril(A(j:m,j), -1) + eye(size(A(j:m,1)));
	      H = (eye(m-j+1,m-j+1) - V * ( lapack_larft( V ) * V' ) );
	      Q(j:m,j:i+n-1) = H * Q(j:m,j:i+n-1);
	end;
*/
