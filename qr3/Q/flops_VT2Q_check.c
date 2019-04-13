#include "qr3.h"

unsigned long int flops_VT2Q_check( int m, int n ){

	unsigned long int flops;

	flops = (( unsigned long int ) 0 );

	flops += flops_ULTinU( n );

	flops += flops_trmm( 'R', m-n, n );

	flops += flops_mLUinA( n );

	flops += (( unsigned long int ) n );

	return flops;

}
