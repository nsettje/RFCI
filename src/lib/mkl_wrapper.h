//Defines C wrappers for FORTRAN MKL functions

#ifndef RFCI_mkl
#define RFCI_mkl
#include <mkl.h>
#include <stddef.h>
#include <math.h>

#define MKL_IP64

//scale a vector x of length n by a constant a
//x->a*x
inline void C_DSCAL(int  n, double a, double *x, int  incx){
	dscal(&n,&a,x,&incx);
}
//multiply two matrices A and B, scale by alpha, and add the result to the scaled matrix beta*C
//C->alpha*A^transa*B^transb+beta*C
inline void C_DGEMM(const char transa, const char transb,int m,int n,int k,double alpha,double* a,int lda, double* b,int ldb, double beta,double* c,int ldc){
	dgemm(&transb,&transa,&n,&m,&k,&alpha,b,&ldb,a,&lda,&beta,c,&ldc);
}
//copy a vector x to a vector y
//y->x
inline void C_DCOPY(int n,const double* x,int incx,double* y,int incy){
	dcopy(&n,x,&incx,y,&incy);
}
//add a scaled vector a*x to a vector y
//y->y+a*x 
inline void C_DAXPY(int n,double a,const double* x,int incx,double* y,int incy){
	daxpy(&n,&a,x,&incx,y,&incy);
}
//take the dot product of two vectors x and y
//return the value
inline double C_DDOT(int n,double* x,int incx,double*y,int incy){
	return ddot(&n,x,&incx,y,&incy);
}
inline void C_DSYEV(const char transa, const char transb,int n, double* a, int lda,double *w, double *work,int lwork){
	int *_info_ = new int;
	dsyev(&transa,&transb,&n,a,&lda,w,work,&lwork,_info_);
	delete _info_;
}
//solves a generalized symmetric-definite eigenproblem
//itype determines the problem to be solved:
//1: A*x = lambda*B*x
//2: A*B*x = lambda*x
//3: B*A*x = lambda*x
inline void C_DSYGVD(int itype,const char jobz,const char uplo,int n,double* a,int lda,double *b,int ldb,double *w, double *work,int lwork,int * iwork, int liwork, int info){
	dsygvd(&itype,&jobz,&uplo,&n,a,&lda,b,&ldb,w,work,&lwork,iwork,&liwork,&info);
}
inline void C_DSYGV(int itype,const char jobz,const char uplo,int n,double* a,int lda,double *b,int ldb,double *w, double *work,int lwork,int * iwork, int liwork, int info){
	dsygv(&itype,&jobz,&uplo,&n,a,&lda,b,&ldb,w,work,&lwork/*,iwork,&liwork*/,&info);
}
inline void C_DGESDD(const char jobz,const int m,int n,double* a,int lda,double *s,double *u,const int ldu,double *vt, const int ldvt,double *work,const int lwork,int * iwork, int info){
	dgesdd(&jobz,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,iwork,&info);
}
#endif
