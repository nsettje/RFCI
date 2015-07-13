//Header file for reorthonormalize_converged_tables.cc
#ifndef reortho_conv_table_hh
#define reortho_conv_table_hh
#include "lib/rfci_generalized_eigen.h"
#include "lib/rfci_wfn.h"
#include "lib/mol_const.h"
#include "lib/wfn.h"
#include "lib/mkl_wrapper.h"
#include <stdlib.h>
#include <string.h>
int reorthonormalize_converged_tables( mol_constant * mol, RFCI_wfn * wfn);
#endif
