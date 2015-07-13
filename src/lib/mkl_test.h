//header file for mkl_test.cc
#ifndef RFCI_mkl_test
#define RFCI_mkl_test
#include "lib/memory.h"
#include "lib/mkl_wrapper.h"
int test_MKL();
int test_C_DDOT();
int test_C_DGEMM();
int test_C_DSYEV();
int test_C_DSYGVD();
#endif
