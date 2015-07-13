#ifndef RFCI_FCI
#define RFCI_FCI
#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>
#include "includes/slater.h" //contains class and functions for explicitly building H and manipulating Slater determinants
#include "includes/permute.h" //contains functions for permuting strings and finding their index in lexical ordering
#endif
void FCI_method(){



/* INITIALIZATIONS */
	//pointers to important constants
	int *nmo_ = new int; //# molecular orbitals
	int *alphae_ = new int; //# alpha electrons
	int *betae_ = new int; //# beta electrons
	char *molname = new char[80], *basisname = new char[80]; //molecule name and basis name

	//find important constants
	//function in MOint.h
	initialize_MO_constants(nmo_,alphae_,betae_,molname,basisname);

	//recast pointers as ints
	int nmo = *(nmo_);
	int alphae = *alphae_;
	int betae = *betae_;

	//initialize one-electron and two-electron integrals
	double **mo_OEI=block_matrix(nmo,nmo); //Pointer to OEI in MO basis
	double *mo_TEI=init_array(nmo*nmo*nmo*nmo);

	//transform MOs, storing in **OEI and *TEI arrays
	MO_transform(mo_OEI,mo_TEI,nmo);

	//Pre-compute hkl' integrals for sigma
	double **hklprime=block_matrix(nmo,nmo); //Pointer to hkl'
	Hklprime(mo_TEI,mo_OEI,hklprime,nmo); //Compute hkl'
	
	//number of alpha strings
	int astringcount = nchoosek(nmo,alphae);
	//number of beta strings
	int bstringcount = nchoosek(nmo,betae);

	int N = astringcount*bstringcount;
/* END INITIALIZATIONS */
	//set to 1 in order to test sigma builds
	int testsigma = 0;
	if(testsigma){
		test_sigma(mo_OEI,hklprime,mo_TEI,alphae, betae, nmo);
	}
	//number of eigenvalues to find
	int roots = 1;
	if(roots>=N){
		roots = N;
	} 
	//Davidson-Liu diagonalization
	printf("Full Dimension = %d\n",astringcount*bstringcount);
	dav_ci_core(roots, hklprime, mo_TEI, mo_OEI, betae, alphae , nmo, molname, basisname,1);
	//Explicit diagonalization of full Hamiltonian
	//printf("\nExact Roots\n");
	/*double **Hfull = block_matrix(N,N);
	build_full_Hamiltonian(Hfull,mo_OEI,mo_TEI,alphae,betae,nmo);
	diagonalize_full_Hamiltonian(Hfull,alphae,betae,nmo,roots);	
	free_block(Hfull);*/
	//Free OEI in MO basis	
	free_block(mo_OEI);

	//free TEI	
	free(mo_TEI);

	//Free hkl'
	free_block(hklprime);
	return Success;
}

//return random double
double randouble(){
	double F = (double)rand() / RAND_MAX;
	return F;
}

void test_excited(int alphae, int nmo){

	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int *Istring = init_int_array(alphae);
	int *Jstring = init_int_array(alphae);
	for(i=0;i<alphae;i++){
		Istring[i]=i;
	}
	for(k=0;k<astringcount;k++){
		for(i=0;i<nmo;i++){
			for(j=0;j<nmo;j++){
				int sgn = excite(Istring,i,j,alphae,nmo,Jstring);
				int indx = stradr(Istring,alphae,nmo);
				for(l=0;l<alphae;l++){
					printf("%d ",Jstring[l]);
				}
				printf("|%d %d|",i,j);
				for(l=0;l<alphae;l++){
					printf("%d ",Istring[l]);
				}
				printf("(%d) [%d]",sgn,indx);
				printf("\n");
			}
		}
		next_combination(Istring,nmo,alphae);
	}
}

void test_sigma(double **mo_OEI,double **hklprime,double *mo_TEI,int alphae, int betae, int nmo){

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N = astringcount*bstringcount;
	double *C = init_array(N);
	double **Hfull = block_matrix(N,N);
	double *sigma = init_array(N);
	double *sigma_vec = init_array(N);
	for(int i=0;i<N;i++){
		C[i] = randouble();
	}
	normalize(C,N);

	get_sigma_vec(C, sigma_vec, mo_TEI, hklprime, alphae, betae, nmo,3);
	build_full_Hamiltonian(Hfull,mo_OEI,mo_TEI,alphae,betae,nmo);
	C_DGEMV('N',N,N,1.0,&(Hfull[0][0]),N,C,1,0.0,sigma,1);
	printf("exact sigma: ");
	for(int i=0;i<N;i++){
		printf("%12.8lf ",sigma[i]);
	}
	printf("\n");
	free_block(Hfull);
	C_DAXPY(N,-1.0,sigma,1,sigma_vec,1);
	printf("Error Norm =%lf\n",C_DDOT(N,sigma_vec,1,sigma_vec,1));
	free(C);
	free(sigma);
	free(sigma_vec);

}

//normalize a vector 
//if the norm is too small, set the vector to zero
double normalize(double *vec, int length){
	double norm = C_DDOT(length,vec,1,vec,1);
	norm = sqrt(norm);
	if(norm>10E-6){
		norm = 1/norm;
		C_DSCAL(length,norm,vec,1);
		return(1/norm);
	}
	else{
		C_DSCAL(length,0,vec,1);
		return 0;
	}
}



// integral terms for sigma, easy to pre-compute
void Hklprime(double *mo_TEI, double **mo_OEI, double **hklprime,int nmo){
	for(int k=0;k<nmo;k++){
		for(int l=0;l<nmo;l++){
			hklprime[k][l]=mo_OEI[k][l];
			for(int j=0;j<nmo;j++){
				hklprime[k][l]-=0.5*mo_TEI[k*nmo*nmo*nmo+j*nmo*nmo+j*nmo+l];

			}
		}
	}
}

//use this one, in general
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


//use this one, in general
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
//takes in C and sigma as vectors
//this works for arbitrary spin states
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

//total sigma build for singlet ground states
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

//Davidson-Liu diagonalization
void dav_ci_core(int M, double **hklprime, double *mo_TEI, double **mo_OEI, int betae, int alphae , int nmo, const char *molname, const char *basisname, int print){

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N = astringcount*bstringcount;

	int max_iters=100;
	int max_dim = M;
	if(8*M<N){
		max_dim = 8*M;
	}
	else{
		max_dim = N;
	}
	int i,j,k,L,I;
	int min_pos, numf, iter, converged, skip_check, init_dim;
	double minimum;
	double BIGNUM = 1E100;
	double cutoff = 10E-9;
	double norm, denom, diff;
	double **sigma=block_matrix(max_dim,N); //sigma vectors, columns
	double **b=block_matrix(max_dim,N); //guess vectors, rows
	double **G=block_matrix(max_dim,max_dim); //approximate Hamiltonian
	double *lambda=init_array(max_dim); //eigenvalues of G
	double *old_lambda=init_array(max_dim); //previous eigenvalues
	double *work=init_array(10*max_dim);

	converged = 0;

	int debug_print = 1;
	
	//build initial C from submatrix guess
	if(N>(max_dim-1)*M){ init_dim = max_dim-M;}
	else if(M<N){
		init_dim = M;
	}
	else{
		M = N;
		init_dim = N;
	}
	printf("Initial Dimension = %d\n",init_dim);
	printf("Maximum Dimension = %d\n",max_dim);
	double **alpha = block_matrix(max_dim,max_dim);
	double **initG = block_matrix(init_dim,init_dim);
		for(i=0; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				initG[i][j] = build_single_Hamiltonian_element(i,j,mo_OEI,mo_TEI,alphae,betae,nmo);
			}
		}
		if(debug_print){
		//fprintf(outfile,"init G P\n");
		//print_mat(initG,init_dim,init_dim,outfile);
		}

		//diagonalize sub-Hamiltonian
		sq_rsp(init_dim, init_dim, initG, lambda, 1, alpha, 1e-12);
		for(i=0;i<init_dim;i++){
			for(j=0;j<init_dim;j++){
				b[i][j] = alpha[j][i];
			}
		}

		
	free_block(initG);
	int built=0;	
	L=init_dim; //current size of subspace
	iter=0;
	converged=0;
//	double **H = block_matrix(N,N);
//	build_full_Hamiltonian(H,mo_OEI,mo_TEI,alphae,betae,nmo);
	while(iter<max_iters && converged < M){
		if(debug_print){
		//fprintf(outfile,"Guess Vectors %d\n",iter);
		//print_mat(b,L,N,outfile);
		}
		if(debug_print){
		double **overlap = block_matrix(L,L);
    		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(b[0][0]), N, 0.0, &(overlap[0][0]), L);
		//fprintf(outfile,"Guess Vector Overlap %d\n",iter); 
		//print_mat(overlap,L,L,outfile);
		free_block(overlap);
		}
		int collapse_check=1;
		printf("\niter = %d\n",iter);
		printf("Subspace Dimension = %d\n",L);
		
	    for(i=built;i<L;i++){
			get_sigma_vec(b[i],sigma[i],mo_TEI,hklprime,alphae,betae,nmo,0);
			//	C_DGEMV('N',N,N,1.0,&(H[0][0]),N,b[i],1,0,sigma[i],1);
}
	    built+=L;
		if(debug_print){
		//fprintf(outfile,"Sigma Vectors %d\n",iter);
		//print_mat(sigma,L,N,outfile);
		}
    	C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(sigma[0][0]), N, 0.0, &(G[0][0]), max_dim); //G = b*sigma
		if(debug_print){
		//fprintf(outfile,"Approximate Hamiltonian %d\n",iter);
		//print_mat(G,L,L,outfile);
		}
	    //diagonalize G
		sq_rsp(L,L, G, lambda, 1, alpha, 1e-12);
		if(debug_print){
		//fprintf(outfile,"Approximate Eigenvectors %d\n",iter);
		//print_mat(alpha,L,L,outfile);
		}
	if(L>=N-1){
	 	converged=M;
		break;
	}	
	double *f; 
	for(k=0;k<M;k++){
	f = init_array(N);
		for(I=0;I<N;I++){
			for(i=0;i<L;i++){
				f[I]+=alpha[k][i]*(sigma[i][I]-lambda[k]*b[i][I]);
			}
		denom=lambda[k]-build_single_Hamiltonian_element(I,I,mo_OEI,mo_TEI,alphae,betae,nmo);
	    		
		if(fabs(denom)>10e-6){ f[I]/=denom;}
    		else{ f[I]=0.0;}
		}
		for(i=0;i<L;i++){
			double proj = C_DDOT(N,f,1,b[i],1);
			C_DAXPY(N,-proj,b[i],1,f,1);
		}
		if(normalize(f,N)<10E6 && L<N){
			for(i=0;i<N;i++){
				b[L][i] = f[i];
			}
			L++;
		}
	free(f);
	}
      	//when L gets too large, collapse subspace
	if(max_dim-L<M){
      		printf("--Collapse--\n");
		//fprintf(outfile,"COLLAPSE\n");
		int new_dim = M;
		if(6*M<max_dim){
			new_dim=6*M;
		}
      		for(i=0;i<L;i++){
      			for(k=0;k<N;k++){
				//STORE OLD GUESS VECTORS IN SIGMA TO SAVE MEMORY
      				sigma[i][k]=b[i][k];//0.0;
				b[i][k] = 0.0;
      			}
      		}
      		for(i=0;i<new_dim;i++){
      			for(j=0;j<L;j++){
      				for(k=0;k<N;k++){
					//THIS SIGMA IS NOT ACTUAL SIGMA
					//It is the old guess vectors stored in sigma to save memory
      					b[i][k]+=alpha[j][i]*sigma[j][k];
      				}
      			}
      		}
		normalize(b[0],N);	
		for(i=1;i<new_dim;i++){	
			normalize(b[i],N);	
			for(j=0;j<i;j++){
				double proj = C_DDOT(N,b[i],1,b[j],1);
				C_DAXPY(N,-proj,b[j],1,b[i],1);
			}
			if(normalize(b[i],N)>10E6){
				for(j=0;j<N;j++){
					b[i][j] = 0.0;
				}
			}
		}
		for(i=0;i<L;i++){
			for(j=0;j<N;j++){
				sigma[i][j] = 0.0;
			}
		}
      		L=new_dim;
      		built=0;
      		collapse_check=0;
      	}
      	//if collapse has not happened, check convergence
      	if(collapse_check){
      		converged = 0;
      		int *conv = init_int_array(M);
      		for(i=0;i<M;i++){
      			conv[i]=0;
      		}
      		for(k=0;k<M;k++){
      			diff=lambda[k]-old_lambda[k];
      			if(fabs(diff)<cutoff){
      				conv[k]=1;
      				converged++;
				      			}
      			printf("%d %2.14lf %2.14lf %1s\n",k,lambda[k],diff,conv[k] == 1 ? "Y" : "N");
      			old_lambda[k]=lambda[k];
      		}
      		free(conv);
      	}
		
	
      	iter++;
      }	
	         
      	//when all roots have converged, calculate final eigenvectors/values
      	
      	if(converged==M){
      		printf("--convergence reached--\ntotal iters = %d\n",iter);
		double *eps=init_array(M);
		double **v=block_matrix(N,M);
		for(i=0;i<M;i++){
			for(k=0;k<N;k++){
				v[k][i]=0.0;
			}
			printf("%d %20.14lf\n",i,lambda[i]);
		}
	
		if(print){
			shared_ptr<Molecule> mol = Process::environment.molecule();
			double rxncoord = 2*mol->fz(3)/1.889726;
			double nuc_rep_energy = mol->nuclear_repulsion_energy();
			char energy_fid[80];
			sprintf(energy_fid,"systems/%s/basis/%s/data/%s_%s_%lf",molname,basisname,molname,basisname,rxncoord);
			printf("Printing to file %s\n",energy_fid);
			FILE *energy_file = fopen(energy_fid,"w");
			fprintf(energy_file,"%20.14lf ",rxncoord);
			fprintf(energy_file,"%20.14lf ",nuc_rep_energy);
			for(i=0;i<M;i++){
				fprintf(energy_file,"%20.14lf ",lambda[i]);
			}
			fprintf(energy_file,"%d\n",iter);
			fclose(energy_file);
		}	
		free_block(v);
		free(eps);
	}		
	      	
      	free_block(b);
      	free_block(sigma);
      	free_block(G);
      	free(old_lambda);
      	//free(lambda);
      	//free(work);
	//free_block(H);
      	printf("Exiting Davidson Routine\n");	
      	
}
}}
