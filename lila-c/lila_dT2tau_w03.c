#include "lila.h"

int lila_dT2tau_w03( int m, int n, int i, int mt, double *T, int ldt, double *tau ){

	double *Tki, *T0j, *Tkk;
	int vb, k, j;

	k   = i % mt;
	Tki = T + k + i*ldt;
	vb  = mt - k; if ( vb > n ) vb = n;

	Tkk = Tki;
	for( k = 0; k < vb; k++){ tau[k] = *Tkk; Tkk=Tkk+1+ldt; }

	j    = i + vb;
	tau += vb;
	T0j  = T + j*ldt;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		Tkk = T0j;
		for( k = 0; k < vb; k++){tau[k] = *Tkk; Tkk=Tkk+1+ldt; }

		j   += vb;
		tau += vb;
		T0j += vb*ldt;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}   
