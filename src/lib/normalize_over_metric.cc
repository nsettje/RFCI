//given a generalized eigenvalue problem AX = BXe, this function computes matrix elements of XBX in order to normalize vectors over the metric
#include "lib/normalize_over_metric.h"

//length = astringcount + 1
double normalize_over_metric_P(int state, double *f, RFCI_wfn *wfn, mol_constant *mol){
	int astringcount = mol->astringcount;
	double *overlap = init_array(astringcount+1);
	get_overlap_P(state,f,overlap,wfn,mol);
	double norm = C_DDOT(astringcount+1,f,1,overlap,1);
	norm = 1/sqrt(norm);
	C_DSCAL(astringcount+1,norm,f,1);
	free(overlap);
	return norm;
}

//length = bstringcount + 1
double normalize_over_metric_Q(int state, double *f, RFCI_wfn *wfn, mol_constant *mol){
	int bstringcount = mol->bstringcount;
	double *overlap = init_array(bstringcount+1);
	get_overlap_Q(state,f,overlap,wfn,mol);
	double norm = C_DDOT(bstringcount+1,f,1,overlap,1);
	norm = 1/sqrt(norm);
	C_DSCAL(bstringcount+1,norm,f,1);
	free(overlap);
	return norm;
}


