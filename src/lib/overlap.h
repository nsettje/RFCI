//Header file for overlap.cc

#ifndef RFCI_overlap
#define RFCI_overlap
#include "lib/mkl_wrapper.h"
#include "lib/permute.h"
#include "lib/rfci_wfn.h"
#include "lib/mol_const.h"
#include "lib/wfn.h"
int get_overlap_P(int state, double *b,double *overlapP,RFCI_wfn *wfn, mol_constant *mol);
int get_overlap_Q(int state, double *b,double *overlapP,RFCI_wfn *wfn, mol_constant *mol);
#endif
