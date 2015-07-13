/* Header file for fci_sigma.cc*/

#ifndef RFCI_fcisigma
#define RFCI_fcisigma
#include "lib/permute.h"
#include "lib/mkl_wrapper.h"
#include "lib/memory.h"
#include <math.h>
#include <stdlib.h>
//#include "util/IO.h"
#include "lib/normalize.h"

//full size FCI sigma build
void get_sigma_vec(double *C, double *sigma, double *mo_TEI, double **hklprime, int alphae, int betae, int nmo,int print);
//supporting functions
void SIGMA1_vec(double *SIGMA, double *Calphabeta, double **hklprime, double *mo_TEI, int alphae, int betae, int nmo );
void SIGMA2_vec(double *SIGMA, double *Calphabeta, double **hklprime, double *mo_TEI, int alphae, int betae, int nmo );
void SIGMA3_vec(double *SIGMA, double *mo_TEI, double *Calphabeta, int alphae, int betae, int astringcount, int bstringcount, int nmo);
#endif
