#include "qr2.h"

int qr2_aux_dlarf_wrapper( char side, int m, int n, double *V, int incv, double tau, double *C, int ldc, double *work){

	dlarf_( &side, &m, &n, V, &incv, &tau, C, &ldc, work);

	return 0;
}
