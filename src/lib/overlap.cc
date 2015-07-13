/* This file contains routines for computing overlap vectors of full length and factored length 
 * for use in computing metric matrices
 *
 */

#include "lib/overlap.h"

int get_overlap_P(int state, double *b,double *overlapP,RFCI_wfn *wfn, mol_constant *mol){
	int debug=0;
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount; 
	int N = astringcount*bstringcount;
	int nterms = wfn->n_terms;
	double *overlap = init_array(astringcount*bstringcount); 
	double *C = init_array(astringcount*bstringcount);
	if(debug){printf("c0 = %lf\n",wfn->c0[state][0]);}
	for(int i=0;i<nterms-1;i++){ 
		expand_factored_wfn(C,1.0/*wfn->c0[state][i]*/,&(wfn->P[state][i*astringcount]),&(wfn->Q[state][i*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	}
	if(debug){
		printf("|psi0> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",C[i]);
		}
		printf("\n");
		printf("| P  > = ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",b[i+1]);
		}
		printf("\n");
		printf("|-Q  > = ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",wfn->Q[state][(nterms-1)*bstringcount+i]);
		}
		printf("\n");
	}
	overlapP[0] = b[0]*C_DDOT(astringcount*bstringcount,C,1,C,1);
	project_expanded_sigmaP(C,&overlapP[1],&(wfn->Q[state][(nterms-1)*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	if(debug){
		printf("P overlap\n");
		for(int i=0;i<astringcount+1;i++){
			printf("%lf ",overlapP[i]);
		}
		printf("\n");
	}
	overlapP[0]+=C_DDOT(astringcount,&overlapP[1],1,&(b[1]),1);
	C_DSCAL(astringcount,b[0],&overlapP[1],1);
	if(debug){
		for(int i=0;i<astringcount+1;i++){
			printf("%lf ",overlapP[i]);
		}
		printf("\n");
	}
	free(C);
	double *D;
	D = init_array(astringcount*bstringcount); 
	if(debug){
		printf("|psiP> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",D[i]);
		}
		printf("\n");
	}
	expand_factored_wfn(D,1.0/*wfn->c0[state][nterms-1]*/,&b[1],&(wfn->Q[state][(nterms-1)*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	if(debug){
		printf("|psiP> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",D[i]);
		}
		printf("\n");
		printf("| P  > = ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",b[i+1]);
		}
		printf("\n");
		printf("|-Q  > = ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",wfn->Q[state][(nterms-1)*bstringcount+i]);
		}
		printf("\n");
	}
	project_expanded_sigmaP(D,&overlapP[1],&(wfn->Q[state][(nterms-1)*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	if(debug){
		for(int i=0;i<astringcount+1;i++){
			printf("%lf ",overlapP[i]);
		}
		printf("\n");
		printf("|c0 P> = ");
		for(int i=0;i<astringcount+1;i++){
			printf("%lf ",b[i]);
		}
		printf("\n");
		printf("P table overlap = %lf\n",C_DDOT(astringcount+1,&overlapP[0],1,&b[0],1));
	}
	free(D);
	free(overlap);
	return 0;
}

int get_overlap_Q(int state, double *b,double *overlapQ,RFCI_wfn *wfn, mol_constant *mol){
	int debug=0;
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount; 
	int N = astringcount*bstringcount;
	int nterms = wfn->n_terms;
	double *overlap = init_array(astringcount*bstringcount); 
	double *C = init_array(astringcount*bstringcount);
	if(debug){printf("c0 = %lf\n",wfn->c0[state][0]);}
	for(int i=0;i<nterms-1;i++){ 
		expand_factored_wfn(C,1.0/*wfn->c0[state][i]*/,&(wfn->P[state][i*astringcount]),&(wfn->Q[state][i*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	}
	if(debug){
		printf("|psi0> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",C[i]);
		}
		printf("\n");
		printf("| P  > = ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",wfn->P[state][(nterms-1)*astringcount+i]);
		}
		printf("\n");
		printf("|-Q  > = ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",b[i+1]);
		}
		printf("\n");
	}
	overlapQ[0] = b[0]*C_DDOT(astringcount*bstringcount,C,1,C,1);
	project_expanded_sigmaQ(C,&overlapQ[1],&(wfn->P[state][(nterms-1)*astringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	if(debug){
		printf("Q overlap\n");
		for(int i=0;i<bstringcount+1;i++){
			printf("%lf ",overlapQ[i]);
		}
		printf("\n");
	}
	overlapQ[0]+=C_DDOT(bstringcount,&overlapQ[1],1,&(b[1]),1);
	C_DSCAL(bstringcount,b[0],&overlapQ[1],1);
	if(debug){
		for(int i=0;i<bstringcount+1;i++){
			printf("%lf ",overlapQ[i]);
		}
		printf("\n");
	}
	free(C);
	double *D;
	D = init_array(astringcount*bstringcount); 
	if(debug){
		printf("|psiQ> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",D[i]);
		}
		printf("\n");
	}
	expand_factored_wfn(D,1.0/*wfn->c0[state][nterms-1]*/,&(wfn->P[state][(nterms-1)*astringcount]),&b[1],*mol->alphae,*mol->betae,*mol->nmo);
	if(debug){
		printf("|psiQ> = ");
		for(int i=0;i<N;i++){
			printf("%lf ",D[i]);
		}
		printf("\n");
		printf("| P  > = ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",wfn->P[state][(nterms-1)*astringcount+i]);
		}
		printf("\n");
		printf("|-Q  > = ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",b[i+1]);
		}
		printf("\n");
	}
	project_expanded_sigmaQ(D,&overlapQ[1],&(wfn->P[state][(nterms-1)*astringcount]),*mol->alphae,*mol->betae,*mol->nmo);
	if(debug){
		printf("|c0 Q> = ");
		for(int i=0;i<bstringcount+1;i++){
			printf("%lf ",b[i]);
		}
		printf("\n");
		printf("Q table overlap = %lf\n",C_DDOT(bstringcount+1,&overlapQ[0],1,&b[0],1));
	}
	free(D);
	free(overlap);
	return 0;
}


