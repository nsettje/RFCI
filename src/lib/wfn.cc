/* This file contains convenient functions for expanding P and Q tables into full size wfns and projectng full size wfns onto P or Q tables
 *
 */

#include "lib/permute.h"

void project_expanded_sigmaP(double *sigma,double *proj_sigma, double *Q, int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int i,j,k,l;
	for(i=0;i<astringcount;i++){
		//printf("P'[%d] = ",i);
		for(j=0;j<bstringcount;j++){
			//printf("C[%d]*Q[%d] +",i*astringcount+j,j);
			proj_sigma[i] += sigma[i*astringcount+j]*Q[j]; 
		}
		//printf("\n");
	}	
	
}

void project_expanded_sigmaQ(double *sigma,double *proj_sigma, double *P, int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int i,j,k,l;
	for(i=0;i<bstringcount;i++){
		//printf("Q'[%d] = ",i);
		for(j=0;j<astringcount;j++){
			//printf("C[%d]*P[%d] +",j*astringcount+i,j);
			proj_sigma[i] += sigma[j*astringcount+i]*P[j]; 
		}
		//printf("\n");
	}	
	
}

void expand_factored_wfn(double *C, double c0, double *P, double *Q, int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int i,j,k,l;
	for(i=0;i<astringcount;i++){
		for(j=0;j<bstringcount;j++){
			C[i*astringcount+j] += c0*P[i]*Q[j];
			//printf("C[%d] = c0*P[%d]*Q[%d] = %lf*%lf*%lf = %lf\n",i*astringcount+j,i,j,c0,P[i],Q[j],C[i*astringcount+j]);
		}
	}	
}

