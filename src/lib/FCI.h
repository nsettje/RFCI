#ifndef RFCI_fci_method_h
#define RFCI_fci_method_h
#include "lib/fci_sigma.h"
#include "lib/mkl_wrapper.h"
#include "lib/hamiltonian.h"
#include "lib/mol_const.h"
#include "lib/memory.h"
#include <stdio.h>
void fci_davidson(mol_constant *mol, char * output, int print);
#endif
