#include "flops.h"

long int flops_geqr3( int int_m, int int_n ){

	long int m, n, flops;

//	flops = (( long int ) 0 );
//	flops += flops_lapack_geqr2( int_m, int_n );
//	flops += flops_larft( int_m, int_n );
//	flops += flops_geqr3_bef_useT( int_n );

	m = ( long int ) int_m;
	n = ( long int ) int_n;

//	flops = (
//
//	         12 * m * n * n    // from GEQR2 
//	       -  4 * n * n * n    // from GEQR2
//	       +  6 * m * n        // from GEQR2
//	       + 34 * n            // from GEQR2
//
//	       +  6 * m * n * n    // from LARFT
//	       -  2 * n * n * n    // from LARFT
//	       -  6 * m * n        // from LARFT
//	       +  3 * n * n        // from LARFT
//	       -  1 * n            // from LARFT
//
//	       +  1 * n * n * n    // from useT
//	       -  3 * n * n        // from useT
//	       +  2 * n            // from useT
//
//	) / 6;

//	flops = (
//
//	         12 * m * n * n    // from GEQR2 
//	       +  6 * m * n * n    // from LARFT
//
//	       -  4 * n * n * n    // from GEQR2
//	       -  2 * n * n * n    // from LARFT
//	       +  1 * n * n * n    // from useT
//
//	       +  6 * m * n        // from GEQR2
//	       -  6 * m * n        // from LARFT
//
//	       +  3 * n * n         // from LARFT
//	       -  3 * n * n         // from useT
//
//	       + 34 * n            // from GEQR2
//	       - 1 * n             // from LARFT
//	       + 2 * n             // from useT
//
//	) / 6;

	flops = (

	         18 * m * n * n

	       -  5 * n * n * n

	       + 35 * n

	) / 6;

	return flops;

}
