/* This file contains functions to initialize RFCI_wfn structs and to read parameters from files
 */

#ifndef RFCI_wfn_memory
#define RFCI_wfn_memory
#include <stdlib.h>
#include "lib/memory.h"
#include "lib/rfci_wfn.h"
#include "lib/mol_const.h"
#include "lib/permute.h"
#include "lib/rfci_settings.h"

struct RFCI_wfn *allocate_memory(rfci_settings *rfci_set, mol_constant *mol){

	int astringcount = mol->astringcount;//nchoosek(nmo,alphae);
	int bstringcount = mol->bstringcount;//nchoosek(nmo,betae);
	int N = astringcount*bstringcount;

	int maximum_davidson_iterations= *rfci_set->RFCI_MAX_DAV_ITERS; 
	int maximum_wavefunction_terms = *rfci_set->RFCI_MAX_WFN_ITERS;
	if(maximum_wavefunction_terms>N){
		maximum_wavefunction_terms = N;
		*rfci_set->RFCI_MAX_WFN_ITERS = N;
	} 
	int roots = *mol->roots;
	struct RFCI_wfn *wfn = (RFCI_wfn*) malloc(sizeof(struct RFCI_wfn));	
	wfn->dav_cutoff= *rfci_set->RFCI_CUTOFF;
	wfn->total_memory_allocated = 0;

	wfn->P = block_matrix(roots,maximum_wavefunction_terms*astringcount);
	wfn->total_memory_allocated += roots*maximum_wavefunction_terms*astringcount*sizeof(double);

	wfn->Q = block_matrix(roots,maximum_wavefunction_terms*bstringcount);
	wfn->total_memory_allocated += roots*maximum_wavefunction_terms*bstringcount*sizeof(double);
//c0 has an extra column because of the indexing scheme for keeping track
//of the values of c0 within a wfn iteration
	wfn->c0 = block_matrix(roots,maximum_wavefunction_terms+1);
	wfn->total_memory_allocated += roots*(maximum_wavefunction_terms+1)*sizeof(double);

	wfn->Ep = init_array(roots);
	wfn->total_memory_allocated += roots*sizeof(int);

	wfn->Eq = init_array(roots);
	wfn->total_memory_allocated += roots*sizeof(int);

	wfn->n_terms = 1; 
	return wfn;
}

void free_rfci_wfn(struct RFCI_wfn *wfn){
	free_block(wfn->P);
	free_block(wfn->Q);
	free_block(wfn->c0);
	free(wfn->Ep);
	free(wfn->Eq);
	wfn->total_memory_allocated = 0; 
}
#endif
