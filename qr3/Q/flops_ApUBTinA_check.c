#include "qr3.h"

unsigned long int flops_ApUBTinA_check( int m, int n ){

	unsigned long int flops;

	int m1, m2, n1, n2;

	flops = (( unsigned long int ) 0 );

	if (( m <= 1 )&&( n<=1 )){

		return (( unsigned long int ) 2 );

	} else {


		m1 = m/2;
		m2 = m-m1;

		n1 = n/2; 
		n2 = n-n1;

		if(( m1 > 0 )&&( n1 > 0 )) flops += flops_ApUBTinA_check( m1, n1 );

		if(( m1 > 0 )&&( n1 > 0 )&&( m2 > 0 )) flops += flops_gemm( m1, n1, m2 );

		if(( m1 > 0 )&&( n2 > 0 )) flops += flops_ApUBTinA_check( m1, n2 );

		if(( m1 > 0 )&&( n2 > 0 )&&( m2 > 0 )) flops += flops_gemm( m1, n2, m2 );

		if(( m2 > 0 )&&( n1 > 0 )) flops += flops_ApUBTinA_check( m2, n1 );

		if(( m2 > 0 )&&( n2 > 0 )) flops += flops_ApUBTinA_check( m2, n2 );

		return flops;

	}

}
