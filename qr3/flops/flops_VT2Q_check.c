#include "flops.h"

long int flops_VT2Q_check( int m, int n ){

	long int flops;

	flops = (( long int ) 0 );

	flops += flops_ULTinU( n );

	flops += flops_trmm( 'R', m-n, n );

	flops += flops_mLUinA( n );

	flops += (( long int ) n );

	return flops;

}
