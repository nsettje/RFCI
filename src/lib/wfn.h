//Header file for wfn.cc
#ifndef RFCI_wavefunction
#define RFCI_wavefunction
void project_expanded_sigmaP(double *sigma,double *proj_sigma, double *Q, int alphae, int betae, int nmo);
void project_expanded_sigmaQ(double *sigma,double *proj_sigma, double *P, int alphae, int betae, int nmo);
void expand_factored_wfn(double *C, double c0, double *P, double *Q, int alphae, int betae, int nmo);
#endif
