//Defines RFCI_wfn struct

#ifndef RFCI_wfn_struct
#define RFCI_wfn_struct
struct RFCI_wfn {
	double **P;
	double **Q;
	double **c0;
	double *Ep;
	double *Eq;
	double *energy;
	int max_dav_iters;
	double dav_cutoff;
	int n_terms;
	int total_memory_allocated;
} ;
#endif
