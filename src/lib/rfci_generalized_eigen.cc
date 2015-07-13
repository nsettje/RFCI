/* This file contains functions necessary for computing all of the components of the generalized eigenproblem.
 *
 * Functions include components to compute:
 * sigma_P
 * sigma_P
 * sigma_P and overlap vector
 * sigma_Q and overlap vector
*/

#include "lib/rfci_generalized_eigen.h"
//#include "lib/wfn.h"
//#include "lib/fci_sigma.h"
void get_sigma_P(int state, double *b,double *sigmaP,RFCI_wfn *wfn, mol_constant *mol){
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount; 
	int N = astringcount*bstringcount;
	double *sigma = init_array(astringcount*bstringcount); 
	double *C = init_array(astringcount*bstringcount); 
	expand_factored_wfn(C,wfn->c0[state][0],&(b[1]),wfn->Q[state],*mol->alphae,*mol->betae,*mol->nmo);
	get_sigma_vec(C,sigma,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
	project_expanded_sigmaP(sigma,&(sigmaP[1]),wfn->Q[state],*mol->alphae,*mol->betae,*mol->nmo);
	free(sigma);
	free(C);
}

void get_sigma_Q(int state, double *b,double *sigmaQ,RFCI_wfn *wfn, mol_constant *mol){
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount; 
	int N = astringcount*bstringcount;
	double *sigma = init_array(astringcount*bstringcount); 
	double *C = init_array(astringcount*bstringcount);
	expand_factored_wfn(C,wfn->c0[state][0],wfn->P[state],&(b[1]),*mol->alphae,*mol->betae,*mol->nmo);
	get_sigma_vec(C,sigma,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
	project_expanded_sigmaQ(sigma,&(sigmaQ[1]),wfn->P[state],*mol->alphae,*mol->betae,*mol->nmo);
	free(sigma);
	free(C);
}

void get_sigma_overlap_P(int state, double *b,double *sigmaP,double *overlapP,RFCI_wfn *wfn, mol_constant *mol){

	int debug = 0;
	
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount; 
	int N = astringcount*bstringcount;
	int nterms = wfn->n_terms;
	double *sigma = init_array(astringcount*bstringcount); 
	double *overlap = init_array(astringcount*bstringcount); 
	double *C = init_array(astringcount*bstringcount);
	//|C> = sum_i c_i |P_iQ_i> for converged tables
	for(int i=0;i<nterms-1;i++){ 
		expand_factored_wfn(C,1.0/*wfn->c0[state][i]*/,&(wfn->P[state][i*astringcount]),&(wfn->Q[state][i*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	}
	if(debug){
		printf("|psi0> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",C[i]);
		}
		printf("\n");
		printf("psi norm = %lf\n",C_DDOT(N,C,1,C,1));
	}

	//sigma = H|C>	
	get_sigma_vec(C,sigma,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
	if(debug){
		printf("|sig0> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",sigma[i]);
		}
		printf("\n");
	}
	//sigmaP[0] = b[0]<C|H|C>
	sigmaP[0] = b[0]*C_DDOT(astringcount*bstringcount,sigma,1,C,1);
	if(debug){
		printf("sigmaP[0] = %lf*%lf = %lf\n",b[0],C_DDOT(astringcount*bstringcount,sigma,1,C,1),sigmaP[0]);
	}
	//sigmaP[1-end] = <C|H|Q>
	project_expanded_sigmaP(sigma,&sigmaP[1],&(wfn->Q[state][(nterms-1)*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	sigmaP[0]+=C_DDOT(astringcount,&sigmaP[1],1,&(b[1]),1);
	if(debug){
		printf("sigmaP[0] = %lf+%lf = %lf\n",sigmaP[0] - C_DDOT(astringcount,&sigmaP[1],1,&(b[1]),1),C_DDOT(astringcount,&sigmaP[1],1,&(b[1]),1),sigmaP[0]);
	}
	C_DSCAL(astringcount,b[0],&sigmaP[1],1);
	free(sigma);
	free(C);
	sigma = init_array(astringcount*bstringcount);
	C = init_array(astringcount*bstringcount);
	// C = |PQ>
	expand_factored_wfn(C,1.0/*wfn->c0[state][nterms-1]*/,&b[1],&(wfn->Q[state][(nterms-1)*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	get_sigma_vec(C,sigma,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
	project_expanded_sigmaP(sigma,&sigmaP[1],&(wfn->Q[state][(nterms-1)*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);

	get_overlap_P(state,b,overlapP,wfn,mol);
	free(C);
	free(sigma);
	free(overlap);
	
}
void get_sigma_overlap_Q(int state, double *b,double *sigmaQ,double *overlapQ,RFCI_wfn *wfn, mol_constant *mol){
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount; 
	int N = astringcount*bstringcount;
	int nterms = wfn->n_terms;
	double *sigma = init_array(astringcount*bstringcount); 
	double *overlap = init_array(astringcount*bstringcount); 
	double *C = init_array(astringcount*bstringcount);
	//|C> = sum_i c_i |P_iQ_i> for converged tables
	for(int i=0;i<nterms-1;i++){ 
		expand_factored_wfn(C,1.0/*wfn->c0[state][i]*/,&(wfn->P[state][i*astringcount]),&(wfn->Q[state][i*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	}
	//sigma = H|C>	
	get_sigma_vec(C,sigma,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
	//sigmaP[0] = b[0]<C|H|C>
	sigmaQ[0] = b[0]*C_DDOT(astringcount*bstringcount,sigma,1,C,1);
	//sigmaP[1-end] = <C|H|Q>
	project_expanded_sigmaQ(sigma,&sigmaQ[1],&(wfn->P[state][(nterms-1)*astringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	sigmaQ[0]+=C_DDOT(bstringcount,&sigmaQ[1],1,&(b[1]),1);
	C_DSCAL(bstringcount,b[0],&sigmaQ[1],1);
	free(sigma);
	free(C);
	sigma = init_array(astringcount*bstringcount);
	C = init_array(astringcount*bstringcount);
	// C = |PQ>
	expand_factored_wfn(C,1.0/*wfn->c0[state][nterms-1]*/,&(wfn->P[state][(nterms-1)*astringcount]),&b[1],*mol->alphae,*mol->betae,*mol->nmo);
	get_sigma_vec(C,sigma,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
	project_expanded_sigmaQ(sigma,&sigmaQ[1],&(wfn->P[state][(nterms-1)*astringcount]),*mol->alphae,*mol->betae,*mol->nmo);

	get_overlap_Q(state,b,overlapQ,wfn,mol);
	free(C);
	free(sigma);
	free(overlap);
	
}


