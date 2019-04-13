#include "qr3.h"

unsigned long int flops_ApUBTinA( int m, int n ){

	unsigned long int flops;

	flops = (( unsigned long int ) m ) * (( unsigned long int ) m ) * (( unsigned long int ) n ) + (( unsigned long int ) m ) * (( unsigned long int ) n );

	return flops;

}
