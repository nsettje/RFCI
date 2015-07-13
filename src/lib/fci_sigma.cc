/*This is the routine and supporting functions for te full-size full CI sigma builds
 *
 * output:
 * double *sigma (values of sigma build, length astringcount*bstringcount)
 *
 * input:
 * double *C (wfn coefficients, length astringcount*bstringcount)
 * mo_OEI (MO one-electron integrals)
 * mo_TEI (MO two-electron integrals)
 * alphae (# alpha electrons)
 * betae  (#  beta electrons)
 * nmo 	  (# MOs)
 * print  (0,1 switch for print to terminal)*/

#include "lib/fci_sigma.h"
//based on the Olsen factoring scheme: 
//sigma_1 - beta-beta
//sigma_2 - alpha-alpha
//sigma_3 - alpha-beta
void get_sigma_vec(double *C, double *sigma, double *mo_TEI, double **hklprime, int alphae, int betae, int nmo,int print){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount=nchoosek(nmo,betae);
	int N = astringcount*bstringcount;
	SIGMA1_vec(sigma, C, hklprime, mo_TEI, alphae, betae, nmo );
	if(print>0){
	printf("sigma 1:     ");
	for(int i=0;i<astringcount*bstringcount;i++){
		printf("%12.8lf ",sigma[i]);
	}
	printf("\n");
	}
	SIGMA2_vec(sigma, C, hklprime, mo_TEI, alphae, betae, nmo );
	if(print>1){
	printf("sigma 1+2:   ");
	for(int i=0;i<astringcount*bstringcount;i++){
		printf("%12.8lf ",sigma[i]);
	}
	printf("\n");
	}
	SIGMA3_vec(sigma,mo_TEI,C,alphae,betae,astringcount,bstringcount,nmo);	
	if(print>2){
	printf("sigma 1+2+3: ");
	for(int i=0;i<astringcount*bstringcount;i++){
		printf("%12.8lf ",sigma[i]);
	}
	printf("\n\n");
	}
		
}

void SIGMA1_vec(double *SIGMA, double *Calphabeta, double **hklprime, double *mo_TEI, int alphae, int betae, int nmo ){
	int kk;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	double * const Fb = new double[bstringcount];	//this is NOT a const double, but a const pointer. Necessary to free Fb.
	int *Ibstring = new int[betae];
	int *Kbstring = new int[betae];
	int *Jbstring = new int[betae];
	for(int i=0;i<betae;i++){
		Ibstring[i]=i;
	}	
	for(int current_string=0;current_string<bstringcount;current_string++){
		for(int i=0;i<bstringcount;i++){
			Fb[i]=0.0;
		}
		for(int k=0;k<nmo;k++){
			for(int l=0;l<betae;l++){
				for(int i=0;i<betae;i++){
					Kbstring[i]=0;
				}
				int sgnkl = excite(Ibstring, k, Ibstring[l], betae, nmo, Kbstring);
				if(sgnkl!=0){
					int Kindex = stradr(Kbstring,betae,nmo);
					Fb[Kindex]+=sgnkl*hklprime[k][Ibstring[l]];
					for(int i=0;i<nmo;i++){
						for(int j=0;j<betae;j++){
							for(int ii=0;ii<betae;ii++){
								Jbstring[ii]=0;
							}
							int sgnij = excite(Kbstring, i, Kbstring[j], betae, nmo, Jbstring); //b
							if(sgnij!=0){
								int Jindex = stradr(Jbstring,betae,nmo);
								Fb[Jindex]+= 0.5*sgnij*sgnkl*mo_TEI[i*nmo*nmo*nmo+Kbstring[j]*nmo*nmo+k*nmo+Ibstring[l]];
							}
						}
					}
				}
			}
		}
			for(int i=0;i<astringcount;i++){
				for(int j=0;j<bstringcount;j++){
					SIGMA[i*bstringcount+current_string]+=Calphabeta[i*bstringcount+j]*Fb[j];
				}
			}	
	
		next_combination(Ibstring,nmo,betae);
		
	}
	delete[] Ibstring;
	delete[] Kbstring;
	delete[] Jbstring;
	delete[] Fb;
}


void SIGMA2_vec(double *SIGMA, double *Calphabeta, double **hklprime, double *mo_TEI, int alphae, int betae, int nmo ){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	double * const Fa = init_array(astringcount);//new double[astringcount]; //this is NOT a const double, but a const pointer. Necessary to free Fa.
	int *const Iastring = new int[alphae];
	int *const Kastring = new int[alphae];
	int *const Jastring = new int[alphae];
	int kk=0;
	int ll=0;
	for(int i=0;i<alphae;i++){
		Iastring[i]=i;
	}	
	for(int current_string=0;current_string<astringcount;current_string++){
		for(int i=0;i<bstringcount;i++){
			Fa[i]=0.0;
		}
		for(int k=0;k<nmo;k++){
			for(int l=0;l<alphae;l++){
				int sgnkl = excite(Iastring, k, Iastring[l], alphae, nmo, Kastring);
				if(sgnkl!=0){
					int Kindex = stradr(Kastring,alphae,nmo);
					Fa[Kindex]+=sgnkl*hklprime[k][Iastring[l]];
					for(int i=0;i<nmo;i++){
						for(int j=0;j<alphae;j++){
							int sgnij = excite(Kastring, i, Kastring[j], alphae, nmo, Jastring); //b
							if(sgnij!=0){
								int Jindex = stradr(Jastring,alphae,nmo);
								Fa[Jindex]+= 0.5*sgnij*sgnkl*mo_TEI[i*nmo*nmo*nmo+Kastring[j]*nmo*nmo+k*nmo+Iastring[l]];
							}
						}
					}
				}
			}
		}

		for(int i=0;i<astringcount;i++){
			for(int j=0;j<bstringcount;j++){
				SIGMA[current_string*bstringcount+j]+=Calphabeta[i*bstringcount+j]*Fa[i];
			}
		}	
	
		
	
					next_combination(Iastring,nmo,alphae);
		
	}
	
	delete[] Iastring;
	delete[] Kastring;
	delete[] Jastring;
  //	delete[] Fa;

}

void SIGMA3_vec(double *SIGMA, double *mo_TEI, double *Calphabeta, int alphae, int betae, int astringcount, int bstringcount, int nmo){
	int *Iastring = init_int_array(alphae);
	int *Ibstring = init_int_array(betae);
	int *Jastring = init_int_array(alphae);
	int *Jbstring = init_int_array(betae);
	for(int i=0;i<alphae;i++){
		Iastring[i]=i;
	}
	
	//loop over Ia
	for(int acount=0;acount<astringcount;acount++){
		//loop over a(Ja)
		for(int k=0;k<nmo;k++){
			for(int l=0;l<alphae;l++){
				int sgnkl = excite(Iastring,k,Iastring[l],alphae,nmo,Jastring);
				if(sgnkl != 0){
					int Jaindex = stradr(Jastring,alphae,nmo);
					for(int i=0;i<betae;i++){
						Ibstring[i]=i;
					}
					//loop over Ib
					for(int bcount=0;bcount<bstringcount;bcount++){
						//loop over b(Jb)
						for(int i=0;i<nmo;i++){
							for(int j=0;j<betae;j++){
								int sgnij = excite(Ibstring,i,Ibstring[j],betae,nmo,Jbstring);
								if(sgnij!=0){
									int Jbindex = stradr(Jbstring,betae,nmo);
									SIGMA[acount*astringcount+bcount]+=sgnij*sgnkl*mo_TEI[i*nmo*nmo*nmo+Ibstring[j]*nmo*nmo+k*nmo+Iastring[l]]*Calphabeta[Jaindex*astringcount+Jbindex];
								}
							}
						} //end loop over Jb
						next_combination(Ibstring,nmo,betae);
					}//end loop over Ib
				}
			}
		}//end loop over Ja
		next_combination(Iastring,nmo,alphae);
	}//end loop over Ia
	
	free(Iastring);
	free(Ibstring);
	free(Jastring);
	free(Jbstring);
}


