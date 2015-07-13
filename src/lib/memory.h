//Header file for memory management functions

#ifndef RFCI_memory
#define RFCI_memory
#include <stdio.h>
double ** block_matrix(int n, int m);
int * init_int_array(int size);
double * init_array(int size);
void free_block(double **array);
void print_mat(double **a, int m, int n, FILE * out);
#endif
