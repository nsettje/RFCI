//Header file for rfci_davidson.cc

#ifndef RFCI_davidson
#define RFCI_davidson
#include "lib/rfci_sigma.h"
#include "lib/normalize.h"
#include "lib/overlap.h"
#include "lib/hamiltonian.h"
#include "lib/mkl_wrapper.h"
#include "lib/rfci_wfn.h"
#include <string.h>
int rfci_davidson(int PorQ, RFCI_wfn *wfn, mol_constant *mol, char * output, int print);
#endif
