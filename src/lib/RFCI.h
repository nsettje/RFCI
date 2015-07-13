//header file for rfci_slow_method.cc

#ifndef RFCI_RFCI_hxx
#define RFCI_RFCI_hxx
#include "lib/mol_const.h"
#include <math.h>
#include "lib/mkl_wrapper.h"
#include "lib/normalize.h"
#include "lib/rfci_davidson.h"
#include "lib/reorthonormalize_converged_tables.h"
int rfci_slow_method(mol_constant *mol, char * output, int print);
#endif
