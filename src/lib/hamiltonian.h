//Header file for buildH.cc

#ifndef RFCI_hamiltonian
#define RFCI_hamiltonian
//#include <util/memory.h>
#include <lib/slater.h>
#include <lib/permute.h>
#include "lib/memory.h"
void build_full_Hamiltonian(double **H, double **mo_OEI, double *mo_TEI, int alphae, int betae, int nmo);
double build_single_Hamiltonian_element(int row, int col, double **mo_OEI, double *mo_TEI, int alphae, int betae, int nmo);
#endif
