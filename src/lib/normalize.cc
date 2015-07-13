
//normalize a vector 
//if the norm is too small, set the vector to zero
#include "lib/normalize.h"
#include <stdio.h>
double normalize(double *vec, int length){
	double norm = C_DDOT(length,vec,1,vec,1);
	norm = sqrt(norm);
	if(norm>10E-9){
		norm = 1/norm;
		C_DSCAL(length,norm,vec,1);
		return(1/norm);
	}
	else{
		C_DSCAL(length,0,vec,1);
		return 0;
	}
}
