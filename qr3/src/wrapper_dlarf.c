#include "qr2.h"

int wrapper_dlarf( char side, int m, int n, double *V, int incv, double tau, double *C, int ldc, double *work){

	dlarf_( &side, &m, &n, V, &incv, &tau, C, &ldc, work);

	return 0;
}
