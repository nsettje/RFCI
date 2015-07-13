//Header file for normalize_over_metric.cc
#ifndef RFCI_normalize_overmetric
#define RFCI_normalize_overmetric
#include "lib/permute.h"
#include "lib/normalize.h"
#include "lib/mkl_wrapper.h"
#include "lib/overlap.h"
#include "lib/rfci_wfn.h"
#include "lib/mol_const.h"
#include <math.h>
//#include "lib/rfci_generalized_eigen.cc"
double normalize_over_metric_P(int state, double *f, RFCI_wfn *wfn, mol_constant *mol);
double normalize_over_metric_Q(int state, double *f, RFCI_wfn *wfn, mol_constant *mol);
#endif
