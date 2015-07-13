//Header file for permute.cc

#ifndef RFCI_permute
#define RFCI_permute
#include <limits.h>
#include "lib/memory.h"
#include <stdio.h>
#include <stdlib.h>
int nchoosek(int n, int k);
void Hklprime(double *mo_TEI, double **mo_OEI, double **hklprime,int nmo);
int excite(int *Istring, int k, int l, int nelec, int nmo, int *excited);     
int excite_fast(int *Istring, int k, int l, int nelec, int nmo, int *excited);
int zindex(int k, int l, int nelec, int nmo);
int stradr(int *str, int nelec, int nmo);
int trunc_stradr(int *str, int nelec, int nmo);
unsigned int next_combination(int * const ar, int n, unsigned int k);
unsigned int next_combination_char(unsigned char * ar, int n,unsigned int k);
int next_combination_fast(int * const ar, int n, int k);
#endif
