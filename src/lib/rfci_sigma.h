//Header for rfci_sigma.cc

#ifndef RFCI_rfcisigma
#define RFCI_rfcisigma
#include "lib/permute.h"
#include "lib/mkl_wrapper.h"
//#include "lib/IO.h"
#include "lib/mol_const.h"
#include "lib/memory.h"
#include <math.h>
int get_NO_augmented_sigmaP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *sigmaP, int n_terms,int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI);
int get_NO_augmented_sigmaQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *sigmaQ, int n_terms,int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI);
int get_factored_sigmaP(double scale, double *P,double *Q,double *Qprime,double *sigmaP, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print);
void get_factored_sigmaQ(double scale, double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print);
void get_factored_sigmaPAA(double *P,double *Q,double *Qprime, double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaPBB(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaPAB(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQAA(double *P,double *Q, double * Pprime, double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQBB(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQAB(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
#endif
