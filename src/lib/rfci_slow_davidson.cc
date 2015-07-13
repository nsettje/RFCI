/* This routine is the workhorse of RFCI. It diagonalizes either P eigenproblem or the Q eigenproblem.
 *
 * output:
 * wfn->P or wfn->Q contains optimized wfn factors
 * wfn->c0 contains optimum normalization parameter
 *
 * input:
 * PorQ (0 = P, 1 = Q)
 * initialized mol_constant struct
 * initialized RFCI_wfn struct
 * output (filename for progress file)
 * print (verbose printing to terminal)*/
#include "lib/rfci_slow_davidson.h"
int rfci_slow_davidson(int PorQ,RFCI_wfn *wfn,mol_constant *mol,char *output,int print){
	//verbose debugging
	int debug_print = 0;
	
	FILE *outfile = fopen(output,"a");

	int i, j, k, I;

	int MAXIT = wfn->max_dav_iters;
	

	int n_terms = wfn->n_terms;

	int astringcount = nchoosek(*mol->nmo,*mol->alphae);
	int bstringcount = nchoosek(*mol->nmo,*mol->betae);
	//set constant parameters
	int M = *mol->roots; //# roots sought
	int N; //dimension of eigenproblem
	char vec_name;
	if(PorQ==1){
		N = astringcount+1;
		vec_name = 'P';	
	}
	else if(PorQ==2){
		N = bstringcount+1;
		vec_name = 'Q';	
	}
	else{
		printf("P or Q parameter must be 1 or 2!\n");
		return 1;
	}
	//maximum subspace dimension
	int maxdim;
	if(8*M<N){
		maxdim = 8*M;
	}
	else{
		maxdim = 2*M;
	}
	if(maxdim < 3){
		maxdim=4;
	}
	//collapse if the subspace gets too big
	int davidson_collapse = 1;

	//speed up convergence for small systems
	if(maxdim<100){
		maxdim = 100;
		davidson_collapse = 0;
	}

	//initial dimension of subspace
	int init_dim;
	if(M > maxdim){
		init_dim = maxdim;
	}
	else if(2*M < maxdim){
		init_dim = 2*M;
	} 
	else{
		init_dim = M;
	}
	init_dim = M;

	//guess vectors		
	double **b = block_matrix(maxdim,N);
	wfn->total_memory_allocated+=maxdim*N*sizeof(double);
	
	//sigma vectors
	double **sigma = block_matrix(maxdim,N);
	wfn->total_memory_allocated+=maxdim*N*sizeof(double);
	
	//overlap vectors for generalized problem (metric)
	double **overlap;
	if(n_terms>1){
		overlap = block_matrix(maxdim,N);
		wfn->total_memory_allocated+=maxdim*N*sizeof(double);
	}
	//copy tables from wfn->P or wfn->Q to b vectors
	if(PorQ==1){
		for(i=0;i<M;i++){
			C_DCOPY(astringcount,&(wfn->P[i][(n_terms-1)*astringcount]),1,&(b[i][1]),1);

			if(n_terms>1){
				b[i][0] = wfn->c0[i][n_terms-2];
			}
		}
	}
	else{
		for(i=0;i<M;i++){
			C_DCOPY(bstringcount,&(wfn->Q[i][(n_terms-1)*bstringcount]),1,&(b[i][1]),1);
			if(n_terms>1){
				b[i][0] = wfn->c0[i][n_terms-2];
			}
		}
				
	}
	int L = init_dim; //current subspace dimension
	int built = 0;//# sigma vectors built so far
	int iter = 0;//iteration counter
	int converged = 0;//# roots converged so far
	int skip_check;//switch to stop convergence check for collapsed iterations
	double *lambda = init_array(maxdim);//energies
	double *old_lambda = init_array(maxdim);//old energies
	double **G=block_matrix(maxdim,maxdim); //approximate Hamiltonian
	double **metric = block_matrix(maxdim,maxdim);//metric matrix
	double **alpha=block_matrix(maxdim,maxdim); //approximate eigenvectors

	if(debug_print){
		fprintf(outfile,"=== DAVIDSON %c %d ===\n",vec_name,wfn->n_terms);
	} 
	while(converged < M && iter < MAXIT){
		//collapse loop is at the top to speed up when many iterations run
		skip_check = 0;
		//COLLAPSE
		if(maxdim-L<M && davidson_collapse){
			if(debug_print){
				fprintf(outfile,"Collapse\n");
				printf("Collapse\n");
			}
			int new_dim = M;
			if(2*M<maxdim){
				//new_dim=2*M;
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
			if(n_terms==1){
				normalize(b[0],N);	
				for(i=1;i<new_dim;i++){	
					normalize(b[i],N);	
					for(j=0;j<i;j++){
						double proj = C_DDOT(N,b[i],1,b[j],1);
						C_DAXPY(N,-proj,b[j],1,b[i],1);
					}
					if(!normalize(b[i],N)){
						for(j=0;j<N;j++){
							b[i][j] = 0.0;
						}
					}
				}
			}
			else{ //reorthonormalize over the metric
				if(PorQ == 1){
					normalize_over_metric_P(0,b[0],wfn,mol);
						if(debug_print){
							double * c_overlap = init_array(N);
							get_overlap_P(0,b[0],c_overlap,wfn,mol);
							print_mat(b,1,N,outfile);
							print_mat(&c_overlap,1,N,outfile);
							free(c_overlap);
						}
					for(i=1;i<new_dim;i++){	
						normalize_over_metric_P(0,b[i],wfn,mol);
						double * collapse_overlap = init_array(N);
						get_overlap_P(0,b[i],collapse_overlap,wfn,mol);
						for(j=0;j<i;j++){
							double proj = C_DDOT(N,collapse_overlap,1,b[j],1);
							C_DAXPY(N,-proj,b[j],1,b[i],1);
						}
						if(!normalize_over_metric_P(0,b[i],wfn,mol)){
							for(j=0;j<N;j++){
								b[i][j] = 0.0;
							}
						}
						if(debug_print){
							free(collapse_overlap);
							collapse_overlap = init_array(N);
							get_overlap_P(0,b[i],collapse_overlap,wfn,mol);
						}
						free(collapse_overlap);
					}
				}
				else{
					normalize_over_metric_Q(0,b[0],wfn,mol);
						if(debug_print){
							double * c_overlap = init_array(N);
							get_overlap_Q(0,b[0],c_overlap,wfn,mol);
							print_mat(b,1,N,outfile);
							print_mat(&c_overlap,1,N,outfile);
							free(c_overlap);
						}
					for(i=1;i<new_dim;i++){	
						normalize_over_metric_Q(0,b[i],wfn,mol);
						double * collapse_overlap = init_array(N);
						get_overlap_Q(0,b[i],collapse_overlap,wfn,mol);
						for(j=0;j<i;j++){
							double proj = C_DDOT(N,collapse_overlap,1,b[j],1);
							C_DAXPY(N,-proj,b[j],1,b[i],1);
						}
						if(!normalize_over_metric_Q(0,b[i],wfn,mol)){
							for(j=0;j<N;j++){
								b[i][j] = 0.0;
							}
						}
						if(debug_print){
							free(collapse_overlap);
							collapse_overlap = init_array(N);
							get_overlap_Q(0,b[i],collapse_overlap,wfn,mol);
						}
						free(collapse_overlap);
					}
				}

			}
			//rezero every array
			for(i=0;i<L;i++){
				for(j=0;j<N;j++){
					sigma[i][j] = 0.0;
					G[i][j] = 0.0;
					if(n_terms>1){
						overlap[i][j] = 0.0;
						alpha[i][j] = 0.0;
					}
				}
			}
			L=new_dim;
			built=0;
			skip_check=1;
		}

		//guess vector overlap matrix
		if(debug_print){
			fprintf(outfile,"b %d\n",iter);
			print_mat(b,L,N,outfile);
			double **boverlap = block_matrix(L,L);
    			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(b[0][0]), N, 0.0, &(boverlap[0][0]), L); 
			fprintf(outfile,"b %c overlap %d\n",vec_name,iter);
			print_mat(boverlap,L,L,outfile);
			free_block(boverlap);
		}

		//BUILD SIGMA 
		//DIAGONALIZE G
		if(n_terms == 1){
			if(PorQ == 1){
				for(i=built;i<L;i++){
					get_sigma_P(0, &(b[i][0]),&(sigma[i][0]),wfn,mol);
					built++;
				}
			}
			else{
				for(i=built;i<L;i++){
					get_sigma_Q(0, &(b[i][0]),&(sigma[i][0]),wfn,mol);
					built++;
				}
			}
			if(debug_print){
				fprintf(outfile,"sigma %c %d\n",vec_name,iter);
				print_mat(sigma,L,N,outfile);
			}
			
    			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim); //G = b*sigma
			sq_rsp(L,L, G, lambda, 1, alpha, 1e-12);
		}
		//GENERALIZED EIGENVALUE PROBLEM
		else{
			if(PorQ == 1){
				for(i=built;i<L;i++){
					get_sigma_overlap_P(0, &(b[i][0]),&(sigma[i][0]),&(overlap[i][0]),wfn,mol);
					built++;
				}
			}
			else{
				for(i=built;i<L;i++){
					get_sigma_overlap_Q(0, &(b[i][0]),&(sigma[i][0]),&(overlap[i][0]),wfn,mol);
					built++;
				}
			}
			if(debug_print){
				fprintf(outfile,"sigma %c %d\n",vec_name,iter);
				print_mat(sigma,L,N,outfile);
				fprintf(outfile,"overlap %c %d\n",vec_name,iter);
				print_mat(overlap,L,N,outfile);
			}
    				C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim); //G = b*sigma	
    				C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N, &(overlap[0][0]), N, 0.0, &(metric[0][0]), maxdim); 
			double *work = init_array(100*L);
			int *iwork = init_int_array(100*L);
			int info = 1;
			fprintf(outfile,"Diagonalizing...\n");
			C_DSYGV(1,'V','L',L,&(G[0][0]),maxdim,&(metric[0][0]),maxdim,lambda,work,100*L,iwork,100*L,info);
			for(i=0;i<L;i++){
				for(j=0;j<L;j++){
					alpha[i][j] = G[j][i];
				}
			}
			fprintf(outfile,"G %c %d\n",vec_name,iter);
			print_mat(G,L,L,outfile);
			fprintf(outfile,"Metric %c %d\n",vec_name,iter);
			print_mat(metric,L,L,outfile);
			fprintf(outfile,"alpha %c %d\n",vec_name,iter);
			print_mat(alpha,L,L,outfile);
			fprintf(outfile,"lambda %c %d\n",vec_name,iter);
			print_mat(&lambda,1,L,outfile);
			free(work);
			free(iwork);
		}

		//TRY TO ADD A NEW  ORTHONORMAL VECTOR
		double *f; 
		for(k=0;k<M;k++){
		f = init_array(N);
			for(I=0;I<N;I++){
				for(i=0;i<L;i++){
					f[I]+=alpha[i][k]*(sigma[i][I]-lambda[k]*b[i][I]);
				}
				double denom=lambda[k]-build_single_Hamiltonian_element(I,I,mol->mo_OEI,mol->mo_TEI,*mol->alphae,*mol->betae,*mol->nmo);
				denom = sqrt(fabs(denom));
				if(fabs(denom)>1e-6){ f[I]/=denom;}
				else{ f[I]=0.0;}
			}
			
			if(n_terms==1){
				for(i=0;i<L;i++){
					double proj = C_DDOT(N,f,1,b[i],1);
					C_DAXPY(N,-proj,b[i],1,f,1);
				}
				if(normalize(f,N)){
					fprintf(outfile,"f %c %d\n",vec_name,iter);
					print_mat(&f,1,N,outfile);	
					for(i=0;i<N;i++){
						b[L][i] = f[i];
					}
					L++;
				}
			}
			else{
				if(PorQ == 1){
					for(i=0;i<L;i++){
						double proj = C_DDOT(N,f,1,overlap[i],1);
						C_DAXPY(N,-proj,b[i],1,f,1);
					}
					if(normalize_over_metric_P(k,f,wfn,mol)<1e6 ){
						fprintf(outfile,"f %c %d\n",vec_name,iter);
						print_mat(&f,1,N,outfile);	
						for(i=0;i<N;i++){
							b[L][i] = f[i];
						}
						L++;
					}
				}
				else{
					for(i=0;i<L;i++){
						double proj = C_DDOT(N,f,1,overlap[i],1);
						C_DAXPY(N,-proj,b[i],1,f,1);
					}
					if(normalize_over_metric_Q(k,f,wfn,mol)<1e6 ){
						fprintf(outfile,"f %c %d\n",vec_name,iter);
						print_mat(&f,1,N,outfile);	
						for(i=0;i<N;i++){
							b[L][i] = f[i];
						}
						L++;
					}

				}
			}
		free(f);
		}

		
		//CHECK CONVERGENCE
		if(!skip_check){
			converged = 0;
			int *conv = init_int_array(M);
			for(i=0;i<M;i++){
				conv[i]=0;
			}
			for(k=0;k<M;k++){
				double diff=lambda[k]-old_lambda[k];
				if(fabs(diff)< wfn->dav_cutoff){
					conv[k]=1;
					converged++;
				}
				old_lambda[k]=lambda[k];
			}
			free(conv);
      		}

		iter++;
	}		
	
	//COMPUTE FINAL EIGENVECTORS
      	if(converged==M){
		if(debug_print){
			printf("%c Converged!----------------------------------\n",vec_name);
			fprintf(outfile,"final alpha\n");
			print_mat(alpha,built,M,outfile);
			fprintf(outfile,"final b\n");
			print_mat(b,built,N,outfile);
		}
		double **v = block_matrix(M,N);
		for(I=0; I < M; I++) {
			for(j=0; j < built; j++) {
				for(i=0; i < N; i++) {
					v[I][i] += alpha[j][I] * b[j][i];
				}
			}
		}
		fprintf(outfile, "final eigenvectors %c\n",vec_name);
		print_mat(v,M,N,outfile);
		char energy_name[200];
		sprintf(energy_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/RFCI/energy/RFCIenergy_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
		FILE *energy_file = fopen(energy_name,"w");
		for(i=0;i<M;i++){
			if(fabs(v[i][0])>1.0){ //if c0 gets too big, renormalize to prevent overflow !!!
				if(n_terms>0){
					if(PorQ==1){	
					normalize(v[i],N);
					}
					else{
					normalize(v[i],N);
					}
				}
				else{
					normalize(v[i],N);
				}	
			}
		}

		//copy tables to wfn struct
		for(i=0;i<M;i++){
			if(PorQ==1){
				wfn->Ep[i]=lambda[i];
				if(n_terms>1){
					wfn->c0[i][n_terms-2] = v[i][0];
				}
				C_DCOPY(astringcount,&(v[i][1]),1,&(wfn->P[i][(n_terms-1)*astringcount]),1);
				if(debug_print){
					fprintf(outfile,"P tables\n");
					print_mat(wfn->P,1,(n_terms+1)*astringcount,outfile);
					fprintf(outfile,"Q tables\n");
					print_mat(wfn->Q,1,(n_terms+1)*bstringcount,outfile);
				}
			}
			else{
				wfn->Eq[i]=lambda[i];
				if(n_terms>1){
					wfn->c0[i][n_terms-2] = v[i][0];
				}
				C_DCOPY(bstringcount,&(v[i][1]),1,&(wfn->Q[i][(n_terms-1)*bstringcount]),1);
				if(debug_print){
					fprintf(outfile,"P tables\n");
					print_mat(wfn->P,1,(n_terms+1)*astringcount,outfile);
					fprintf(outfile,"Q tables\n");
					print_mat(wfn->Q,1,(n_terms+1)*bstringcount,outfile);
				}
			}
		}
		

		free_block(v);
	}
	else{
		printf("DAVIDSON DID NOT CONVERGE!\n"); //print to terminal if routine fails
	}
	if(n_terms>1){
		free_block(overlap);
		wfn->total_memory_allocated-=maxdim*N*sizeof(double);
	}
	
			free(metric);
	free_block(b);
	wfn->total_memory_allocated-=maxdim*N*sizeof(double);
	
	free_block(sigma);
	wfn->total_memory_allocated-=maxdim*N*sizeof(double);

	fclose(outfile);

	return 0;
} 
