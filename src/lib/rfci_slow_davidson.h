#ifndef RFCI_rfci_slow_dav
#define RFCI_rfci_slow_dav
#include "lib/rfci_generalized_eigen.h"
#include "lib/sq_rsp.h"
#include "lib/hamiltonian.h"
#include "lib/normalize_over_metric.h"
int rfci_slow_davidson(int PorQ,RFCI_wfn *wfn,mol_constant *mol,char *output,int print);
#endif
