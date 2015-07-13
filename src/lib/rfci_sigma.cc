/* This file contains factored sigma builds for P and Q tables
 */

#include "lib/rfci_sigma.h"


int get_NO_augmented_sigmaP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *sigmaP, int n_terms,int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI){
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	//loop over converged P and Q tables
	double *sigbuff;
	sigmaP[0] = 0.0;
	for(i=0;i<n_terms-1;i++){
		//compute <psi0|H|psi0>
		for(j=0;j<n_terms-1;j++){ //loop over converged P and Q tables
			sigbuff = init_array(astringcount);
			//sigma(P_i,Q_i,Q_j)
			get_factored_sigmaP(1.0,&P[i*astringcount],&Q[i*astringcount],&Q[j*astringcount],sigbuff,alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			//c_i*c_j*sigma(P_i,Q_i,Q_j)*P_j
			sigmaP[0]+=Pnew[0]*c0[i]*c0[j]*C_DDOT(astringcount,sigbuff,1,&P[j*astringcount],1);
			free(sigbuff);
		}
		//compute <psi0|H|aQ>
		sigbuff = init_array(astringcount);
		get_factored_sigmaP(1.0,&P[i*astringcount],&Q[i*astringcount],Qnew,sigbuff,alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		//update sigmaP[0]
		sigmaP[0]+=c0[i]*C_DDOT(astringcount,sigbuff,1,&Pnew[1],1);
		//copy scaled sigbuff to sigmaP
		C_DAXPY(astringcount,c0[i],sigbuff,1,&sigmaP[1],1);
		free(sigbuff);
	}
	//multiply sigmaP by c0
	C_DSCAL(astringcount,Pnew[0],&sigmaP[1],1);
	//compute <a'Q|H|PQ>
	get_factored_sigmaP(1.0,&Pnew[1],Qnew,Qnew,&sigmaP[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	return 0;
}

int get_NO_augmented_sigmaQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *sigmaQ, int n_terms,int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI){

	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	//loop over converged P and Q tables
	double *sigbuff;
	sigmaQ[0] = 0.0;
	for(i=0;i<n_terms-1;i++){
		//compute <psi0|H|psi0>
		for(j=0;j<n_terms-1;j++){ //loop over converged P and Q tables
			sigbuff = init_array(bstringcount);
			//sigma(P_i,Q_i,Q_j)
			get_factored_sigmaQ(1.0,&P[i*astringcount],&Q[i*bstringcount],&P[j*astringcount],sigbuff,alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			//c_i*c_j*sigma(P_i,Q_i,Q_j)*P_j
			sigmaQ[0]+=Qnew[0]*c0[i]*c0[j]*C_DDOT(bstringcount,sigbuff,1,&Q[j*bstringcount],1);
			free(sigbuff);
		}
		//compute <psi0|H|aQ>
		sigbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&P[i*astringcount],&Q[i*astringcount],Pnew,sigbuff,alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		//update sigmaP[0]
		sigmaQ[0]+=c0[i]*C_DDOT(bstringcount,sigbuff,1,&Qnew[1],1);
		//copy scaled sigbuff to sigmaP
		C_DAXPY(bstringcount,c0[i],sigbuff,1,&sigmaQ[1],1);
		free(sigbuff);
	}
	//multiply sigmaP by c0
	C_DSCAL(bstringcount,Qnew[0],&sigmaQ[1],1);
	//compute <a'Q|H|PQ>
	get_factored_sigmaQ(1.0,Pnew,&Qnew[1],Pnew,&sigmaQ[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	return 0;

}



//compute sigma^P for a given state up to a given number of alpha and beta terms in the trial wavefunction, multiply each element of the final vector by the factor "scale"
int get_factored_sigmaP(double scale, double *P,double *Q,double *Qprime,double *sigmaP, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print){
	int astringcount = nchoosek(nmo,alphae);
	print = 3;
	get_factored_sigmaPAA(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>0){
		printf("sigmaP 1:     ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
	get_factored_sigmaPBB(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>1){
		printf("sigmaP 1+2:   ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}

	get_factored_sigmaPAB(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>2){
		printf("sigmaP 1+2+3: ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
	if(fabs(scale-1.0)>10E-6){
		C_DSCAL(astringcount,scale,sigmaP,1);
		printf("scaledsigmaP: ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
	return 1;	
}



//compute sigma^Q for a given state up to a given number of alpha and beta terms in the trial wavefunction
void get_factored_sigmaQ(double scale, double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print){
	int bstringcount = nchoosek(nmo,betae);
	get_factored_sigmaQBB(P,Q,Pprime,sigmaQ, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>1){
		printf("sigmaQ 2:     ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}

	get_factored_sigmaQAA(P,Q,Pprime,sigmaQ, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>0){
		printf("sigmaQ 1+2:   ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}
	get_factored_sigmaQAB(P,Q,Pprime,sigmaQ, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>2){
		printf("sigmaQ 1+2+3: ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}
	if(fabs(scale-1.0)>10E-6){
		C_DSCAL(bstringcount,scale,sigmaQ,1);
		printf("scaledsigmaQ: ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}
}

//compute alpha-alpha component of sigma^P for given state
void get_factored_sigmaPAA(double *P,double *Q,double *Qprime, double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;
	if(debug_print){
		printf("PAA:\n");
	}
	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Iaprimestring = init_int_array(alphae);
	int *Iastring = init_int_array(alphae);
	int *Kastring = init_int_array(alphae);
	int *Jastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	double QB_prefactor = 0.0;
	//pre-compute sum_beta Q_beta*Q_beta
	QB_prefactor+=C_DDOT(bstringcount,Qprime,1,Q,1);
	if(fabs(QB_prefactor)>10E-6){	
		//loop up to highest excitation included in alpha wfn
		for(L=0;L<astringcount;L++){
			//loop over alpha 
			for(int kk=0;kk<alphae;kk++){
				Iastring[kk]=kk;
			}
			for(I=0;I<astringcount;I++){
				//first excitations
				for(l=0;l<alphae;l++){
					//one-electron coupling
					//only excite up to strings that can match Iaprime
					for(k=0;k<alphae;k++){
						int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
						if(sgnkl != 0 ){
							int Kindex = stradr(Kastring,alphae,nmo);
							if(Kindex == L){
								sigmaP[L]+=sgnkl*mo_OEIprime[Iaprimestring[k]][Iastring[l]]*P[I]*QB_prefactor;
							}
						}
					}
					//two-electron coupling
					for(k=0;k<nmo;k++){
						int sgnkl = excite(Iastring,k,Iastring[l],alphae,nmo,Kastring);
						if(sgnkl != 0){
							//second excitations
							for(i=0;i<alphae;i++){
								for(j=0;j<nmo;j++){
									int sgnij = excite(Kastring,Iaprimestring[i],j,alphae,nmo,Jastring);
									if(sgnij != 0 ){
										int Jindex = stradr(Jastring,alphae,nmo);
										if(Jindex == L){
											sigmaP[L]+=0.5*QB_prefactor*sgnkl*sgnij*mo_TEI[Iaprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Iastring[l]]*P[I];
										}
									}
								}
							}
						}
					}
				}
				next_combination(Iastring,nmo,alphae);	
			}
			next_combination(Iaprimestring,nmo,alphae);
		}
		if(debug_print){
			printf("\n");
		}
	}
	free(Iaprimestring);
	free(Iastring);
	free(Kastring);
	free(Jastring);
}

//compute beta-beta component of sigma^P
void get_factored_sigmaPBB(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	

	int debug_print = 0;
	if(debug_print){
		printf("PBB:\n");
	}

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibstring = init_int_array(betae);
	int *Ibprimestring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Jbstring = init_int_array(betae); 

	double Qfactor = 0.0;	

	for(i=0;i<betae;i++){
		Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	//loop over beta prime
	for(I=0;I<bstringcount;I++){
		//loop over beta
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		}
		for(J=0;J<bstringcount;J++){
			//first excitations
			for(l=0;l<betae;l++){
				//one-electron coupling
				for(k=0;k<betae;k++){
					int sgnkl = excite(Ibstring,Ibprimestring[k],Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						int Kindex = stradr(Kbstring,betae,nmo);
						if(Kindex == I){
							Qfactor+=sgnkl*mo_OEIprime[Ibprimestring[k]][Ibstring[l]]*Qprime[I]*Q[J];

						}
					}
				}
				//two-electron coupling
				for(k=0;k<nmo;k++){
					int sgnkl = excite(Ibstring,k,Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						//second excitations
						for(i=0;i<betae;i++){
							for(j=0;j<nmo;j++){
								int sgnij = excite(Kbstring,Ibprimestring[i],j,betae,nmo,Jbstring);
								int Jindex = stradr(Jbstring,betae,nmo);
								if(Jindex == I && sgnij != 0){
									Qfactor+=0.5*sgnkl*sgnij*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Ibstring[l]]*Qprime[I]*Q[J];
								}
							}
						}
					}
				}
			}

			next_combination(Ibstring,nmo,betae);	
		}
		next_combination(Ibprimestring,nmo,betae);	
	}
	for(i=0;i<astringcount;i++){
		sigmaP[i]+=P[i]*Qfactor;
	}
	if(debug_print){
		printf("\n");
	}
	free(Ibprimestring);
	free(Ibstring);
	free(Kbstring);
	free(Jbstring);
}

//compute alpha-beta component of sigma^P
void get_factored_sigmaPAB(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){


	int debug_print = 0;
	if(debug_print){
		printf("PAB:\n");
	}

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibprimestring = init_int_array(betae);
	int *Ibstring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Iastring = init_int_array(alphae); 
	int *Iaprimestring = init_int_array(alphae); 
	int *Kastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	for(i=0;i<betae;i++){
		Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	//initialize beta ij density matrix
	double **beta_density = block_matrix(nmo,nmo);
	//loop over beta prime
	for(I=0;I<bstringcount;I++){	
		//loop over beta
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		} 
		for(J=0;J<bstringcount;J++){
			//excitation ij
			for(i=0;i<betae;i++){
				for(j=0;j<betae;j++){
					//only loop over i in the bra string and j in the ket string
					int sgnij = excite(Ibstring,Ibprimestring[i],Ibstring[j],betae,nmo,Kbstring);
					int Kindex = stradr(Kbstring,betae,nmo);
					if(Kindex == I && sgnij != 0){
						beta_density[Ibprimestring[i]][Ibstring[j]]+=sgnij*Qprime[I]*Q[J];
					}
				}
			}
			next_combination(Ibstring,nmo,betae);
		}
		next_combination(Ibprimestring,nmo,betae);
	}
	//print_mat(beta_density,nmo,nmo,outfile);
	//highest excitation for alpha'
	for(J=0;J<astringcount;J++){
		//loop over alpha
		for(int kk=0;kk<alphae;kk++){
			Iastring[kk]=kk;
		} 
		for(I=0;I<astringcount;I++){
			//excitations kl
			for(k=0;k<alphae;k++){
				for(l=0;l<alphae;l++){
					int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
					if(sgnkl != 0 ){
						int Kindex = stradr(Kastring,alphae,nmo);
						if(Kindex == J ){
							for(i=0;i<nmo;i++){
								for(j=0;j<nmo;j++){
									sigmaP[J]+=beta_density[i][j]*sgnkl*mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+Iaprimestring[k]*nmo+Iastring[l]]*P[I];
								}
							}
						}
					}
				}
			}
			next_combination(Iastring,nmo,alphae);
		}
		next_combination(Iaprimestring,nmo,alphae);
	}
	if(debug_print){
		printf("\n");
	}
	free_block(beta_density);
	free(Ibprimestring);
	free(Ibstring);
	free(Kbstring);
	free(Iastring);
	free(Iaprimestring);
	free(Kastring);
}

//compute alpha-alpha component of sigma^Q
void get_factored_sigmaQAA(double *P,double *Q, double * Pprime, double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Iastring = init_int_array(alphae);
	int *Iaprimestring = init_int_array(alphae);
	int *Kastring = init_int_array(alphae);
	int *Jastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	//loop over alpha prime
	for(I=0;I<astringcount;I++){
		//loop over alpha
		for(int kk=0;kk<alphae;kk++){
			Iastring[kk]=kk;
		}
		for(J=0;J<astringcount;J++){
			//first excitations
			for(l=0;l<alphae;l++){
				//one-electron coupling
				for(k=0;k<alphae;k++){
					int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
					if(sgnkl != 0 ){
						int Kindex = stradr(Kastring,alphae,nmo);
						if(Kindex == I){
							for(K=0;K<bstringcount;K++){
								sigmaQ[K]+=sgnkl*mo_OEIprime[Iaprimestring[k]][Iastring[l]]*Pprime[I]*P[J]*Q[K];
							}
						}
					}
				}
				//two-electron coupling
				for(k=0;k<nmo;k++){
					int sgnkl = excite(Iastring,k,Iastring[l],alphae,nmo,Kastring);
					if(sgnkl != 0){
						//second excitations
						for(i=0;i<alphae;i++){
							for(j=0;j<nmo;j++){
								int sgnij = excite(Kastring,Iaprimestring[i],j,alphae,nmo,Jastring);
								if(sgnij != 0){
									int Jindex = stradr(Jastring,alphae,nmo);
									if(Jindex == I){
										for(K=0;K<bstringcount;K++){
											sigmaQ[K]+=0.5*sgnkl*sgnij*mo_TEI[Iaprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Iastring[l]]*Pprime[I]*P[J]*Q[K];
										}
									}
								}
							}
						}
					}
				}
			}
			next_combination(Iastring,nmo,alphae);	
		}
		next_combination(Iaprimestring,nmo,alphae);	
	}
	if(debug_print){
		printf("\n");
	}
	free(Iastring);
	free(Iaprimestring);
	free(Kastring);
	free(Jastring);
}


//compute beta-beta component of sigma^Q
void get_factored_sigmaQBB(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibprimestring = init_int_array(betae);
	int *Ibstring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Jbstring = init_int_array(betae); 
	for(i=0;i<betae;i++){
		//Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	double PA_prefactor = 0.0;
	//pre-compute sum_beta Q_beta*Q_beta
	PA_prefactor=C_DDOT(astringcount,Pprime,1,P,1);
	for(L=0;L<bstringcount;L++){
		//loop over beta
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		}
		for(I=0;I<bstringcount;I++){
			//first excitations
			for(l=0;l<betae;l++){
				//one-electron coupling
				for(k=0;k<betae;k++){
					int sgnkl = excite(Ibstring,Ibprimestring[k],Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						int Kindex = stradr(Kbstring,betae,nmo);
						if(Kindex == L){
							sigmaQ[L]+=PA_prefactor*sgnkl*mo_OEIprime[Ibprimestring[k]][Ibstring[l]]*Q[I];
						}								
					}
				}
				//two-electron coupling
				for(k=0;k<nmo;k++){
					int sgnkl = excite(Ibstring,k,Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						//second excitations
						for(i=0;i<betae;i++){
							for(j=0;j<nmo;j++){
								int sgnij = excite(Kbstring,Ibprimestring[i],j,betae,nmo,Jbstring);
								if(sgnij != 0){
									int Jindex = stradr(Jbstring,betae,nmo);
									if(Jindex == L){
										sigmaQ[L]+=PA_prefactor*0.5*sgnkl*sgnij*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Ibstring[l]]*Q[I];
																										}
								}
							}
						}
					}
				}
			}
			next_combination(Ibstring,nmo,betae);	
		}
		next_combination(Ibprimestring,nmo,betae);
	}
	if(debug_print){
		printf("\n");
	}
	free(Ibstring);
	free(Ibprimestring);
	free(Kbstring);
	free(Jbstring);
}

//compute alpha-beta component of sigma^Q
void get_factored_sigmaQAB(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibprimestring = init_int_array(betae);
	int *Ibstring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Iastring = init_int_array(alphae); 
	int *Iaprimestring = init_int_array(alphae); 
	int *Kastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	for(i=0;i<betae;i++){
		Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	//loop over alpha prime	
	//initialize alpha ij density matrix
	double **alpha_density = block_matrix(nmo,nmo);
	//loop over alpha prime
	for(I=0;I<astringcount;I++){	
		//loop over alpha 
		for(int kk=0;kk<alphae;kk++){
			Iastring[kk]=kk;
		}
		for(J=0;J<astringcount;J++){
			//excitation ij
			for(i=0;i<alphae;i++){
				for(j=0;j<alphae;j++){
					//only loop over i in the bra string and j in the ket string
					int sgnij = excite(Iastring,Iaprimestring[i],Iastring[j],alphae,nmo,Kastring);
					int Kindex = stradr(Kastring,alphae,nmo);
					if(Kindex == I && sgnij != 0){
						alpha_density[Iaprimestring[i]][Iastring[j]]+=sgnij*Pprime[I]*P[J];
					}
				}
			}
			next_combination(Iastring,nmo,alphae);
		}
		next_combination(Iaprimestring,nmo,alphae);
	}
	//highest excitation for beta'
	for(L=0;L<bstringcount;L++){
		//loop over beta 
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		}
		for(I=0;I<bstringcount;I++){
			//excitations kl
			for(k=0;k<betae;k++){
				for(l=0;l<betae;l++){
					int sgnkl = excite(Ibstring,Ibprimestring[k],Ibstring[l],betae,nmo,Kbstring);
					int Kindex = stradr(Kbstring,betae,nmo);
					if(Kindex == L && sgnkl != 0 ){
						for(i=0;i<nmo;i++){
							for(j=0;j<nmo;j++){
								sigmaQ[L]+=alpha_density[i][j]*sgnkl*mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+Ibprimestring[k]*nmo+Ibstring[l]]*Q[I];
								if(debug_print){
								//	printf("+ %1.1f*%d*(%d %d|%d %d)Q[%d] ",alpha_density[i][j],sgnkl,i,j,Ibprimestring[k],Ibstring[l],I);
								}
							}
						}
					}
				}
			}
			next_combination(Ibstring,nmo,betae);
		}
		next_combination(Ibprimestring,nmo,betae);
	}
	if(debug_print){
		printf("\n");
	}
	free_block(alpha_density);


	free(Ibstring);
	free(Ibprimestring);
	free(Kbstring);
	free(Iastring);
	free(Iaprimestring);
	free(Kastring);
}


