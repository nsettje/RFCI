#include <iostream>
#include <mkl.h>
#include <stddef.h>
#include <math.h>
#include "/home/settje/FCI/lib/block_matrix.cc"
#include "/home/settje/FCI/lib/init_array.cc"
#include "/home/settje/FCI/lib/print_mat.cc"

void print_test();

static inline void C_DSCAL(int  n, const double a, double *x, int  incx); 

int test()
{
	using namespace std;
	print_test();
	double ** A = block_matrix(10,10,false);
	double ** B = block_matrix(10,10,false);
	double ** C = block_matrix(10,10,false);
	A[0][0] = 1.0;
	B[0][0] = 2.0;
	double * x = init_array(10);
	double * y = init_array(10);
	//dscal();
	x[0]=1.0;
	C_DSCAL(10,2.0,x,1);
	printf("x0 = %lf\n",x[0]);
	free(x);
	free(y);
	free_block(A);
	free_block(B);
	free_block(C);
	return 0;
}

void print_test(){
	std::cout << "Hello World\n";
}

static inline void C_DSCAL(int  n, const double a, double *x, int  incx){
	dscal(&n,&a,x,&incx);
}

