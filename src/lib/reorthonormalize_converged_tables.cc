/* This function takes a set of P and Q tables and diagonalizes with respect to the overlap metric. 
 * This returns linear combinations of original tables that corresponds to an eigenbasis.
 *
 */

#include "lib/reorthonormalize_converged_tables.h"
int reorthonormalize_converged_tables( mol_constant * mol, RFCI_wfn * wfn){

	int i,j,k,l;
	int state = 0;

	int n_terms = wfn->n_terms;
	//n_terms++;
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount;

	double **sigma = block_matrix(n_terms,n_terms);
	double **metric= block_matrix(n_terms,n_terms);

	double *Ci = init_array(astringcount*bstringcount);
	double *Cj = init_array(astringcount*bstringcount);
	double *SIGMA = init_array(astringcount*bstringcount);
	
	//COMPUTE SIGMA AND METRIC MATRICES
	for(i=0;i<n_terms;i++){ 
		expand_factored_wfn(Ci,1.0,&(wfn->P[state][i*astringcount]),&(wfn->Q[state][i*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
		get_sigma_vec(Ci,SIGMA,mol->mo_TEI,mol->mo_OEIprime,*mol->alphae,*mol->betae,*mol->nmo,0);
		for(j=0;j<n_terms;j++){
			expand_factored_wfn(Cj,1.0,&(wfn->P[state][j*astringcount]),&(wfn->Q[state][j*bstringcount]),*mol->alphae,*mol->betae,*mol->nmo);
			sigma[i][j] = C_DDOT(astringcount*bstringcount,Cj,1,SIGMA,1);
			metric[i][j] = C_DDOT(astringcount*bstringcount,Ci,1,Cj,1);
			memset(Cj,0.0,sizeof(double)*astringcount*bstringcount);
		}
		memset(Ci,0.0,sizeof(double)*astringcount*bstringcount);
		memset(SIGMA,0.0,sizeof(double)*astringcount*bstringcount);
	}
	//PRINT MATRICES
	/*printf("SIGMA\n");
	for(i=0;i<n_terms;i++){
		for(j=0;j<n_terms;j++){
			printf("%lf ",sigma[i][j]);
		}
		printf("\n");
	}
	printf("METRIC\n");
	for(i=0;i<n_terms;i++){
		for(j=0;j<n_terms;j++){
			printf("%lf ",metric[i][j]);
		}
		printf("\n");
	}*/
	//DIAGONALIZE
	double *work = init_array(100*n_terms);
	int *iwork = init_int_array(100*n_terms);
	double *lambda = init_array(n_terms);
	int info = 1;
	C_DSYGV(1,'V','L',n_terms,&(sigma[0][0]),n_terms,&(metric[0][0]),n_terms,lambda,work,100*n_terms,iwork,100*n_terms,info);
	for(i=0;i<n_terms;i++){
		normalize(&(sigma[i][0]),n_terms);
	}
	free(work);
	free(iwork);
	//PRINT MATRICES
	/*printf("SIGMA\n");
	for(i=0;i<n_terms;i++){
		for(j=0;j<n_terms;j++){
			printf("%lf ",sigma[i][j]);
		}
		printf("\n");
	}
	printf("METRIC\n");
	for(i=0;i<n_terms;i++){
		for(j=0;j<n_terms;j++){
			printf("%lf ",metric[i][j]);
		}
		printf("\n");
	}*/
	
	double *Pbuff = init_array(astringcount);
	double *Qbuff = init_array(bstringcount);
	//APPLY EIGENVECTORS TO P,Q TABLES
	for(i=0;i<n_terms;i++){
		for(j=0;j<n_terms;j++){
			C_DAXPY(astringcount,sigma[i][j],&wfn->P[state][j*astringcount],1,Pbuff,1);
			C_DAXPY(bstringcount,sigma[i][j],&wfn->Q[state][j*bstringcount],1,Qbuff,1);
		}
		C_DCOPY(astringcount,Pbuff,1,&wfn->P[state][i*astringcount],1);
		C_DCOPY(bstringcount,Qbuff,1,&wfn->Q[state][i*bstringcount],1);
		memset(Pbuff,0.0,sizeof(double)*astringcount);
		memset(Qbuff,0.0,sizeof(double)*bstringcount);
	}


	free(Pbuff);
	free(Qbuff);
	free_block(sigma);	
	free_block(metric);
	free(Ci);	
	free(Cj);	
	free(SIGMA);	
	return 0;
}

