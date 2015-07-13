//Header file for normalize.cc and normalize_over_metric.cc
#ifndef RFCI_normalize
#define RFCI_normalize
#include "lib/mkl_wrapper.h"
double normalize(double *vec, int length);
double normalize_over_metricP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, int n_terms,int alphae, int betae, int nmo);
double normalize_over_metricQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, int n_terms,int alphae, int betae, int nmo);
#endif
