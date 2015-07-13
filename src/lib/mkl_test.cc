//File with small tests for each MKL function

#include "lib/mkl_test.h"

int test_MKL(){
	int errors = 0;
	errors+=test_C_DDOT();
	errors+=test_C_DGEMM();
	errors+=test_C_DSYEV();
	errors+=test_C_DSYGVD();
	return errors;
}

int test_C_DSYGVD(){
	int dim = 2;
	int error = 0;
	double **X = block_matrix(dim,dim);
	double **A = block_matrix(dim,dim);
	for(int i=0;i<dim;i++){
		A[i][i] = 1.0;
		for(int j=0;j<dim;j++){
			X[i][j] = i + j + 1.0;
		}
	}
	X[0][0] = -0.3173662;
	X[1][0] = 0.7248762;
	X[0][1] = X[1][0];
	X[1][1] = -1.6556441;

	A[0][0] = 0.1916874;
	A[1][0] = 0.4378212;
	A[0][1] = A[1][0];
	A[1][1] = 1.000000;

	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			printf("%lf ",X[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			printf("%lf ",A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	double *lambda = init_array(dim);
	double *work = init_array(50*dim);
	int *iwork = init_int_array(20*dim);
	int info = 0;
	C_DSYGV(1,'V','U',dim,&(X[0][0]),dim,&(A[0][0]),dim,lambda,work,50*dim,iwork,20*dim,info);
	printf("info = %d\n",info);
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			printf("%lf ",X[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			printf("%lf ",A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	free(work);
	free(lambda);
	free(iwork);
	free_block(X);
	free_block(A);
	return 0;
}

int test_C_DGEMM(){
	int error = 0;
	double **X = block_matrix(3,3);
	double **A = block_matrix(3,3);
	A[0][0] = 1.0;
	A[1][1] = 1.0;
	A[2][2] = 1.0;
	C_DGEMM('N','T',3,3,3,1.0,&(A[0][0]),3,&(A[0][0]),3,0.0,&(X[0][0]),3);
	if(3.0!=A[0][0]+A[1][1]+A[2][2]){
		error = 1;
	}
	free_block(X);
	free_block(A);
	return error;

}

int test_C_DSYEV(){
	int error = 0;
	double **A = block_matrix(3,3);
	A[0][0] = 1.0;
	A[1][1] = 1.0;
	A[2][2] = 1.0;
	double *lambda = init_array(3);
	double *work = init_array(100);
	C_DSYEV('V','U',3,A[0],3,lambda,work,100);
	free(work);
	free(lambda);
	free_block(A);
	return error;
}

int test_C_DDOT(){
	int error = 0;
	double *X = init_array(2);
	X[0] = -1;
	X[1] = 1;
	if(2.0!=C_DDOT(2,X,1,X,1)){
		error = 1;
	}
	free(X);
	return error;
}
