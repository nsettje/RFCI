//Header file for rfci_generalized_eigen.cc

#ifndef RFCI_rfci_gen_eigen
#define RFCI_rfci_gen_eigen
#include "lib/mol_const.h"
#include "lib/rfci_wfn.h"
#include "lib/fci_sigma.h"
#include "lib/memory.h"
#include "lib/permute.h" 
#include "lib/wfn.h" 
#include "lib/overlap.h" 
void get_sigma_overlap_P(int state, double *b,double *sigmaP,double *overlapP,RFCI_wfn *wfn, mol_constant *mol);
void get_sigma_overlap_Q(int state, double *b,double *sigmaP,double *overlapP,RFCI_wfn *wfn, mol_constant *mol);
void get_sigma_P(int state, double *b,double *sigmaP,RFCI_wfn *wfn, mol_constant *mol);
void get_sigma_Q(int state, double *b,double *sigmaQ,RFCI_wfn *wfn, mol_constant *mol);
#endif
