/* This is the diagonalization routine for full CI. It is a Davidson iterative subspace diagonalization.
 *
 * output: 
 * prints converged energies to a file (flat format)
 * prints wfn coefficients to a different file (flat format)
 *
 * input:
 * mol_constant object that contains various electronic constants
 * char * output (file name for progress file)
 * print (0 or 1 to switch off or on verbose printing)*/

#include "lib/FCI.h"
#include "lib/sq_rsp.h"
#include "lib/fci_settings.h" //macros defining all variables beginning with "FCI_"
#ifndef FCI_DEBUG
#define FCI_DEBUG
#endif



void fci_davidson(mol_constant *mol, char * output, int print){
	int i,j,k,I;
	//struct fci_settings;	
	int M = *mol->roots;
	struct fci_settings *settings = read_fci_settings();
	//output file
	FILE * out = fopen(output,"w");

	//collect constants from mol
	int nmo = *mol->nmo;
	int alphae= *mol->alphae;
	int betae = *mol->betae;
	
	//compute relevant 
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N = astringcount*bstringcount;

	//print the full dimension to terminal
	printf("N: %d\n",N);

	//if you have fewer than 1000 determinants, don't bother to collapse
	int FCI_DAVIDSON_COLLAPSE = *settings->FCI_DAVIDSON_COLLAPSE;

	//maximum number of iterations
	int max_iters = *settings->FCI_MAX_ITERS;
	
	double cutoff = *settings->FCI_CUTOFF;
	free_fci_settings(settings);
	double norm, denom, diff;

	//verbose printing
	int debug_print = 1;
	
	//maximum dimension of subspace
	int max_dim = M;
	if(16*M<N){
		max_dim =16*M;
	}
	else{
		max_dim = N;
	}
	//initial dimension of subspace	
	int init_dim;
	if(N>(max_dim-1)*M){ init_dim = max_dim-M;}
	else if(M<N){
		init_dim = M;
	}
	else{
		M = N;
		init_dim = N;
	}
	if(8*M<max_dim){init_dim = 8*M;}
	else{init_dim = M;}
	
	printf("Initial Dimension = %d\n",init_dim);
	printf("Maximum Dimension = %d\n",max_dim);

	if(max_dim == N){
		FCI_DAVIDSON_COLLAPSE = 0;
	}
	
	//build sub-Hamiltonian to diagonalize for an initial guess
	double **G=block_matrix(max_dim,max_dim); //approximate Hamiltonian
		for(i=0; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				G[i][j] = build_single_Hamiltonian_element(i,j,mol->mo_OEI,mol->mo_TEI,alphae,betae,nmo);
			}
		}
		fprintf(out,"G init\n");
		print_mat(G,init_dim,init_dim,out);
	//diagonalize sub-Hamiltonian
	double **b=block_matrix(max_dim,N); //guess vectors, rows
	double *lambda=init_array(max_dim); //eigenvalues of G
	double **alpha = block_matrix(max_dim,max_dim);
	sq_rsp(init_dim, init_dim, G, lambda, 1, alpha, 1e-12);

	for(i=0;i<init_dim;i++){
		for(j=0;j<init_dim;j++){
			b[i][j] = alpha[j][i]; //transpose
		}
	}
		
	int built = 0;	//vectors built this iteration
	int L = init_dim; //current size of subspace
	int iter = 0; //iteration index
	int converged = 0; //number of converged eigenstates
	int skip_check; //flag to skip convergence check if subspace collapses this iteration
	double **sigma=block_matrix(max_dim,N); //sigma vectors
	double *old_lambda=init_array(max_dim); //previous eigenvalues
	
//BEGIN DAVIDSON ITERATIONS
	while(iter<max_iters && converged < M){
		
		printf("\niter = %d\n",iter);
		printf("Subspace Dimension = %d\n",L);
		
		int collapse_check=1;

		//for debugging
		if(debug_print){

			fprintf(out,"Guess Vectors %d\n",iter);
			print_mat(b,L,N,out);

			double **overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(b[0][0]), N, 0.0, &(overlap[0][0]), L);
			fprintf(out,"Guess Vector Overlap %d\n",iter); 
			print_mat(overlap,L,L,out);
			free_block(overlap);
		}

		//build sigma vectors
		for(i=built;i<L;i++){
			printf("Building sigma %d\n",i);
			get_sigma_vec(b[i],sigma[i],mol->mo_TEI,mol->mo_OEIprime,alphae,betae,nmo,0);
		}
	    	built+=L;
		if(debug_print){
			fprintf(out,"Sigma Vectors %d\n",iter);
			print_mat(sigma,L,N,out);
		}
		//build sub-Hamiltonian
    		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(sigma[0][0]), N, 0.0, &(G[0][0]), max_dim); //G = b*sigma
		fprintf(out,"G %d\n",iter);
		print_mat(G,L,L,out);
		//diagonalize sub-Hamiltonian
		sq_rsp(L,L, G, lambda, 1, alpha, 1e-12);
		//try to add new vectors to the space
		double *f; 
		for(k=0;k<M;k++){
		f = init_array(N);
			for(I=0;I<N;I++){
				for(i=0;i<L;i++){
					f[I]+=alpha[k][i]*(sigma[i][I]-lambda[k]*b[i][I]);
				}
				denom=lambda[k]-build_single_Hamiltonian_element(I,I,mol->mo_OEI,mol->mo_TEI,alphae,betae,nmo);
				if(fabs(denom)>10e-6){ f[I]/=denom;}
				else{ f[I]=0.0;}
			}
			for(i=0;i<L;i++){
				double proj = C_DDOT(N,f,1,b[i],1);
				C_DAXPY(N,-proj,b[i],1,f,1);
			}
			if(normalize(f,N)<10E6 && L<max_dim){
				for(i=0;i<N;i++){
					b[L][i] = f[i];
				}
				L++;
			}
		free(f);
		}

		//when L gets too large, if the collapse flag is set in fci_settings.txt then collapse subspace
		if(max_dim-L<M && FCI_DAVIDSON_COLLAPSE){
			printf("--Collapse--\n");
			//printf("Collapse if not add = %d\n",collapse_if_not_add); 
			if(debug_print){fprintf(out,"COLLAPSE\n");}
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
				if(normalize(b[i],N)<10E-6){
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
		char rfci_energy_name[200];
		sprintf(rfci_energy_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/FCI/energy/FCIenergy_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
		FILE *rfci_energy_file = fopen(rfci_energy_name,"w+");
		if(rfci_energy_file!=NULL){
			printf("Writing energies to: %s\n",rfci_energy_name);
			for(i=0;i<M;i++){
				fprintf(rfci_energy_file,"%20.16lf\n",lambda[i]);
			}
		}
		if(rfci_energy_file!=NULL){
		fclose(rfci_energy_file);
		}

		double **v=block_matrix(N,M);
		for(i=0;i<M;i++){
			for(j=0;j<L;j++){
				for(k=0;k<N;k++){
					v[k][i]+=alpha[j][i]*b[j][k];
				}
			}
		}
				char wfn_name[200];
		sprintf(wfn_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/FCI/wfn/FCIwfn_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
		FILE *wfn_file = fopen(wfn_name,"w");
		if(wfn_file!=NULL){
				printf("Printing wavefunction coefficients to: %s\n",wfn_name);
				for(i=0;i<N;i++){
					for(j=0;j<M;j++){
						fprintf(wfn_file,"%20.14lf ",v[i][j]);
					}
					fprintf(wfn_file,"\n");
				}
							//}
		}
		fclose(wfn_file);
		free_block(v);
	}		
	      	
	fclose(out);
      	free_block(b);
      	free_block(sigma);
      	free_block(G);
      	free(old_lambda);
      	free(lambda);
      	printf("Exiting Davidson Routine\n");	
      	
}


