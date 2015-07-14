/* This routine solves the RFCI eigenproblem using the factored sigma builds
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

#include "lib/rfci_davidson.h"
int rfci_davidson(int PorQ, RFCI_wfn *wfn, mol_constant *mol, char * output, int print){
	int betae = *mol->betae;
	int alphae= *mol->alphae;
	int nmo = *mol->nmo;
	int n_terms = wfn->n_terms;
	int i, j, k, L, I;
	int iter, *conv, converged, maxdim, skip_check, built;
	int init_dim;
	double **b,  **sigma, **G;
	double *lambda, **alpha, *f, *old_lambda,**eig_overlap;
	double norm, denom, diff;
	double BIGNUM = 10E100;
	int MAXIT = 1000;
	char vec_name;
	FILE *outfile = fopen(output,"a");
	double cutoff = 10E-7;
	double **P = wfn->P;
	double **Q = wfn->Q;
	double **c0 = wfn->c0;
	double **mo_OEI = mol->mo_OEI;
	double **mo_OEIprime = mol->mo_OEIprime;
	double *mo_TEI = mol->mo_TEI;

	double *energy = wfn->energy;

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N;

	int M = 1;
 1;
	//set sigma vector length for P or Q
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

	if(M>N){
		printf("Requesting more roots than exist! Solving for maximum number of roots.\n");
		M = N;
	}

	//set maximum dimension for subspace
	if(16*M<N){
		maxdim = 16*M;
	}
	else{
		maxdim = 2*M;
	}
	if(maxdim < 3){
		maxdim=3;
	}
	printf("maxdim = %d\n",maxdim);
	b = block_matrix(maxdim, N); //guess vectors, stored by row 
	sigma = block_matrix(maxdim,N); // sigma vectors, stored by row
	G = block_matrix(maxdim, maxdim); // Davidson mini-Hamitonian
	alpha = block_matrix(maxdim, maxdim); // eigenvectors of G 
	lambda = init_array(maxdim); //eigenvalues of G
	old_lambda = init_array(maxdim); // approximate roots from previous iteration 

	//set initial dimension for subspace	
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
	//copy current best guesses to Davidson guesses
	//scale guesses by normalization factors
	if(PorQ==1){ //P guesses
		for(i=0;i<M;i++){
			C_DCOPY(astringcount,&(P[i][(n_terms-1)*astringcount]),1,&(b[i][1]),1);
			if(n_terms>1){
				b[0][0] = c0[i][n_terms-1];
			}
		}
	}
	else{ //Q guesses
		for(i=0;i<M;i++){
			C_DCOPY(bstringcount,&(Q[i][(n_terms-1)*bstringcount]),1,&(b[i][1]),1);
			if(n_terms>1){
				b[0][0] = c0[i][n_terms-1];
			}
		}
	}

	//if starting with more guess vectors than best guesses, use sub-Hamiltonian eigenvectors as remaining Davidson guesses
	if(init_dim>M){
		printf("Using sub-Hamiltonian eigenvectors as remaining initial guesses\n");
		//build sub-Hamiltonian element-by-element
		double **initG = block_matrix(init_dim,init_dim);
		if(PorQ==1){ //P sub-Hamiltonian
			for(i=0; i < init_dim; i++) {
				for(j=0; j < init_dim; j++){
					initG[i][j] = build_single_Hamiltonian_element((n_terms-1)*astringcount+i,(n_terms-1)*astringcount+j,mo_OEI,mo_TEI,alphae,betae,nmo);
				}
			}
		}
		else{ //Q sub-Hamiltonian
			for(i=0; i < init_dim; i++) {
				for(j=0; j < init_dim; j++){
					initG[i][j] = build_single_Hamiltonian_element(i*astringcount+n_terms-1,j*astringcount+n_terms-1,mo_OEI,mo_TEI,alphae,betae,nmo);
				}
			}
		}
		
		//diagonalize sub-Hamiltonian
		double *work = init_array(10*init_dim);
		double *eig= init_array(10*init_dim);
		//C_DSYEV('N','U',init_dim,init_dim,1,eig,work,10*init_dim);
		free(work);
		free(eig);

		//copy eigenvectors to Davidson guesses
		for(i=M;i<init_dim;i++){
			for(j=0;j<init_dim;j++){
				b[i][j+1] = alpha[j][i];
			}
		}
		//orthonormalize new guess vectors to previous guess vectors
		for(i=M;i<init_dim;i++){
			for(j=0;j<i;j++){
				double proj = C_DDOT(N,b[i],1,b[j],1);
				C_DAXPY(N,-proj,b[j],1,b[i],1);
			}
			normalize(b[i],N);
		}

	}

	//set current dimension L to initial dimension
	L = init_dim;

	//initialize iterations and number of converged roots
	iter =0;
	converged = 0;
	conv = init_int_array(M); 

	//keep track of current number of sigma vectors in memory
	built = 0;

	if(n_terms>1){
		eig_overlap = block_matrix(maxdim,N);
	}

	//FOR DEBUGGING
	if(PorQ == 2){
		//fprintf(outfile,"P table\n");
		//print_mat(&P[0],1,astringcount,outfile);
	}
	//END DEBUGGING

	//Build subspace and diagonalize until all roots converge or reach maximum iterations
	while(converged < M && iter < MAXIT) {
		printf("Iteration %d\n",iter);
		//keep track of whether this iteration required subspace collapse
		skip_check = 0;

		//print guess vectors (b are guesses)
		//fprintf(outfile,"B %d (L = %d)\n",iter,L);
		//print_mat(b,L,N,outfile);
	
		//check orthonormality of b vectors
		int print_overlap = 1;
		if(print_overlap){
			double **overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(b[0][0]), N, 0.0, &(overlap[0][0]), L);
			//fprintf(outfile,"%c overlap\n",vec_name);
			print_mat(overlap,L,L,outfile);
			free_block(overlap);	
		}
		//build sigma vectors and overlap vectors
		if(n_terms>1){
			if(PorQ==1){
				for(i=built;i<L;i++){
					get_NO_augmented_sigmaP(&(c0[0][0]),&(b[i][0]), &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), &(sigma[i][0]), n_terms,alphae, betae, nmo, mo_OEIprime, mo_TEI);
					//get_overlapP(&(c0[0][0]),&(b[i][0]), &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
					printf("Built sigma %d\n",built);
					built++;
				}
			}
			else{
				for(i=built;i<L;i++){
					get_NO_augmented_sigmaQ(&(c0[0][0]),&(P[0][(n_terms-1)*astringcount]),&(b[i][0]), &(P[0][0]), &(Q[0][0]), &(sigma[i][0]), n_terms,alphae, betae, nmo, mo_OEIprime, mo_TEI);
					//get_overlapP(&(c0[0][0]), &(b[i][0]),&(P[0][(n_terms-1)*astringcount]),&(Q[0][0]),&(P[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
					printf("Built sigma %d\n",built);
					built++;
				}
			}
		}
		else{
			if(PorQ==1){
				for(i=built;i<L;i++){
					get_factored_sigmaP(1.0,&(b[i][1]),&(Q[0][0]),&(Q[0][0]),&(sigma[i][1]), alphae, betae, nmo, mo_OEIprime,  mo_TEI,0);
					printf("Built sigma %d\n",built);
					built++;
				}
			}
			else{
				for(i=built;i<L;i++){
					get_factored_sigmaQ(1.0,&(P[0][0]),&(b[i][1]),&(P[0][0]),&(sigma[i][1]), alphae, betae, nmo, mo_OEIprime,  mo_TEI,0);
					printf("Built sigma %d\n",built);
					built++;
				}
			}
		}
		printf("Built Sigma\n");
		//print sigma vectors
		//fprintf(outfile,"sigma %c %d\n",vec_name,iter);
		//print_mat(sigma,built,N,outfile);
	
		/* form mini-matrix */
		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim);


		//print quadratic forms
		//fprintf(outfile,"subH%c %d\n",vec_name,iter);
		//print_mat(G,L,L,outfile);
	
	
		if(n_terms>1){
			//fprintf(outfile,"%c Overlap %d\n",vec_name,iter);
			//print_mat(eig_overlap,L,N,outfile);
			//print_mat(b,L,N,outfile);
			double **G_overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(eig_overlap[0][0]), N, 0.0, &(G_overlap[0][0]), L);
			double **B_overlap = block_matrix(L,L);
			C_DCOPY(L*L,&(G_overlap[0][0]),1,&(B_overlap[0][0]),1);
			
			//fprintf(outfile,"%c Generalized Overlap %d\n",vec_name,iter);
			//print_mat(G_overlap,L,L,outfile);
			double *work = init_array(50*L);
			int *iwork = init_int_array(20*L);
			//fprintf(outfile,"B matrix before Diagonalization\n");
			//print_mat(G_overlap,L,L,outfile);
			int info = 0;
			C_DSYGVD(1,'V','U',L,&(G[0][0]),maxdim,&(G_overlap[0][0]),L,lambda,work,50*L,iwork,20*L,info);
			//fprintf(outfile,"B matrix after Diagonalization\n");
			//print_mat(G_overlap,L,L,outfile);
			//double **eigenvector_overlap = block_matrix(L,L);
			//C_DGEMM('n','t', L, L, L, 1.0, &(G[0][0]), maxdim,&(G[0][0]), maxdim, 0.0, &(eigenvector_overlap[0][0]), L);
			//fprintf(outfile,"%c Eigenvector Overlap %d\n",vec_name,iter);
			//print_mat(eigenvector_overlap,L,L,outfile);
			//free_block(eigenvector_overlap);
			//C_DSYEV('V','U',L,&(G[0][0]),maxdim,lambda,work,10*astringcount*bstringcount);
			for(i=0;i<L;i++){
				for(j=0;j<L;j++){
					alpha[i][j] = G[j][i];
				}
			}
			free(work);
			free(iwork);
			double **temp = block_matrix(L,L);
			print_mat(B_overlap,L,L,outfile);
			C_DGEMM('n','n', L, L, L, 1.0, &(B_overlap[0][0]), L,&(alpha[0][0]), maxdim, 0.0, &(temp[0][0]), L);
			
			C_DGEMM('t','n', L, L, L, 1.0, &(alpha[0][0]), maxdim,&(temp[0][0]), L, 0.0, &(G_overlap[0][0]), L);
			//fprintf(outfile,"Normalization through Metric\n");
			//print_mat(G_overlap,L,L,outfile);
			free_block(temp);
			free_block(G_overlap);
		}
		else{
			/* diagonalize mini-matrix */
			//sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);
		double *Lwork = init_array(10*init_dim);
		double *Leig= init_array(10*init_dim);
		printf("Minimatrix: %f\n",G[0][0]);
		//C_DSYEV('N','U',L,L,1,lambda,,10*l);
		free(Lwork);
		free(Leig);
		}

		//fprintf(outfile,"%c eigenvalues\n",vec_name);
		for(j=0;j<M;j++){
			//fprintf(outfile,"%lf\n",lambda[j]);
		}
		//fprintf(outfile,"subH%c eigenvectors\n",vec_name);
		//print_mat(alpha,L,L,outfile);

		/* form preconditioned residue vectors */
		for(k=0; k < M; k++) {//rows
			f = init_array(N);
			for(I=0; I < N; I++) { //cols
				for(i=0; i < L; i++) {
					f[I] += alpha[i][k] * (sigma[i][I] - lambda[k] * b[i][I]);
				}
				denom=lambda[k]-build_single_Hamiltonian_element(I,I,mo_OEI,mo_TEI,alphae,betae,nmo);
				if(fabs(denom) > 10e-6) {
					f[I] /= denom;
				}
				else{
					f[I] = 0.0;
				}
			}
			//orthonormalize vectors wrt current subspace
			if(n_terms == 1 ){
				for(I=0;I<L;I++){
					double proj = C_DDOT(N,b[I],1,f,1);
					C_DAXPY(N,-proj,b[I],1,f,1);
				}
			}
			else{
				for(I=0;I<L;I++){
					double proj = C_DDOT(N,eig_overlap[I],1,f,1);
					C_DAXPY(N,-proj,b[I],1,f,1);
				}
			}
			//if normalization constant of new vectors greater than tolerance, increase size of subspace
			if(n_terms == 1){
				//if(L<maxdim){
					normalize(f,N);
					C_DCOPY(N,f,1,b[L],1);
					L++;
				//}
			}
			else{
				if(L<maxdim){
					if(PorQ == 1){
						if(10E-6 < normalize_over_metricQ(&(c0[0][0]),&(Q[0][(n_terms-1)*bstringcount]), f,&(Q[0][0]), &(P[0][0]), n_terms,alphae, betae, nmo)){
							C_DCOPY(N,f,1,b[L],1);
							L++;
						}
					}
					else{
						if(10E-6 < normalize_over_metricQ(&(c0[0][0]),&(P[0][(n_terms-1)*astringcount]), f,&(P[0][0]), &(Q[0][0]), n_terms,alphae, betae, nmo)){
							C_DCOPY(N,f,1,b[L],1);
							L++;
						}
					}
				}
			}
			
			free(f);
		}

	
	//when L gets too large, collapse subspace
	if(maxdim-L<M){
      		printf("--Collapse--\n");
		//fprintf(outfile,"COLLAPSE\n");
		int new_dim = M;
		printf("M = %d\n",N);
		if(4*M<maxdim){
			new_dim=4*M;
		}
      		for(i=0;i<L;i++){
      			for(k=0;k<N;k++){
				//STORE OLD GUESS VECTORS IN SIGMA TO SAVE MEMORY
      				sigma[i][k]=b[i][k];//0.0;
				b[i][k] = 0.0;
      			}
      		}
      		for(i=0;i<new_dim;i++){
      			for(j=0;j<built;j++){
      				for(k=0;k<N;k++){
					//THIS SIGMA IS NOT ACTUAL SIGMA
					//It is the old guess vectors stored in sigma to save memory
      					b[i][k]+=alpha[j][i]*sigma[j][k];
      				}
      			}
		}
      		for(i=0;i<new_dim;i++){
			//normalize new vectors over the metric
			if(n_terms == 1){
				normalize(b[0],N);	
			}
			else{
				if(PorQ == 1){
					normalize_over_metricQ(&(c0[0][0]),&(Q[0][(n_terms-1)*bstringcount]), b[i],&(Q[0][0]), &(P[0][0]), n_terms,alphae, betae, nmo);
				}
				else{
					normalize_over_metricQ(&(c0[0][0]),&(P[0][(n_terms-1)*astringcount]), b[i],&(P[0][0]), &(Q[0][0]), n_terms,alphae, betae, nmo);
				}
			}

		}
		if(n_terms == 1){
			for(i=1;i<new_dim;i++){	
				normalize(b[i],N);	
				for(j=0;j<i;j++){
					double proj = C_DDOT(N,b[i],1,b[j],1);
					C_DAXPY(N,-proj,b[j],1,b[i],1);
				}
				normalize(b[i],N);	
				/*if(normalize(b[i],N)>10E6){
					for(j=0;j<N;j++){
						b[i][j] = 0.0;
					}
				}*/
			}
		}
		else{
			if(PorQ == 1){
				for(i=1;i<new_dim;i++){	
					//get_overlapP(&(c0[0][0]),&(b[i][0]), &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
					for(j=0;j<i;j++){
						double proj = C_DDOT(N,eig_overlap[i],1,b[j],1);
						C_DAXPY(N,-proj,b[j],1,b[i],1);
					}
					normalize_over_metricQ(&(c0[0][0]),&(Q[0][(n_terms-1)*bstringcount]), b[i],&(Q[0][0]), &(P[0][0]), n_terms,alphae, betae, nmo);
				}
			}
			else{
				for(i=1;i<new_dim;i++){	
					//get_overlapP(&(c0[0][0]), &(b[i][0]),&(P[0][(n_terms-1)*astringcount]),&(Q[0][0]),&(P[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
					for(j=0;j<i;j++){
						double proj = C_DDOT(N,eig_overlap[i],1,b[j],1);
						C_DAXPY(N,-proj,b[j],1,b[i],1);
					}
					normalize_over_metricQ(&(c0[0][0]),&(P[0][(n_terms-1)*astringcount]), b[i],&(P[0][0]), &(Q[0][0]), n_terms,alphae, betae, nmo);
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
		skip_check = 1;
      	}
      	//if collapse has not happened, check convergence
      	if(!skip_check){
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
		
	
		if(PorQ==2){
			//return 0;	
		}	
      	iter++;
      }	
	         
      	//when all roots have converged, calculate final eigenvectors/values
      	
      	if(converged==M){
      		printf("--convergence reached--\ntotal iters = %d\n",iter);
		for(i=0;i<M;i++){
			//energy[i] = lambda[i];
			printf("%d %20.14lf\n",i,lambda[i]);
		}
	printf("HERE\n");
		//fprintf(outfile,"Final B\n");
		//print_mat(b,L,N,outfile);
		//fprintf(outfile,"Final Alpha\n");
		//print_mat(alpha,L,M,outfile);
		double final_L = L;
		if(n_terms == 1){final_L = L-1;}
		double **v=block_matrix(M,N);
		for(I=0; I < M; I++) {
			for(j=0; j < final_L; j++) {
				for(i=0; i < N; i++) {
					v[I][i] += alpha[j][I] * b[j][i];
				}
			}
		}
		//fprintf(outfile,"%c Eigenvectors\n",vec_name);
		//print_mat(v,M,N,outfile);
		if(PorQ==1){
			for(i=0;i<M;i++){
				if(n_terms == 1){
					//normalize(v[0],N);
				}
				else{
					//normalize_over_metricP(&(c0[0][0]),v[0], &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), n_terms,alphae, betae, nmo);
				}
				if(n_terms>1){c0[i][n_terms-1] = v[i][0];}
				C_DCOPY(astringcount,&(v[i][1]),1,&(P[i][(n_terms-1)*astringcount]),1);
			}
		}
		else{
			for(i=0;i<M;i++){
				normalize_over_metricQ(&(c0[0][0]),&(P[0][(n_terms-1)*astringcount]), v[0],&(P[0][0]), &(Q[0][0]), n_terms,alphae, betae, nmo);
				if(n_terms>1){c0[i][n_terms-1] = v[i][0];}
				C_DCOPY(astringcount,&(v[i][1]),1,&(Q[i][(n_terms-1)*bstringcount]),1);
			}
		}
		//fprintf(outfile,"Normalization Factors\n");
		//print_mat(c0,1,n_terms,outfile);
		//fprintf(outfile,"Eigenvector Norm = %lf\n",C_DDOT(N,&(P[0][0]),1,&(P[0][0]),1));
		//FOR DEBUGGING
		/*if(PorQ == 1 && n_terms>1){
			double *v_overlap = init_array(astringcount +1);
			get_overlapP(&(c0[0][0]),&(v[0][0]), &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), v_overlap, n_terms,alphae, betae, nmo);
			double final_overlap = C_DDOT(astringcount+1,v[0],1,v_overlap,1);
			fprintf(outfile,"Final Eigenvector Overlap = %lf\n",C_DDOT(astringcount+1,v[0],1,v_overlap,1));
			if(fabs(1.0-final_overlap)>10E-6){
				printf("Wavefunction is not normalized: norm = %lf\n",final_overlap);
			}
			free(v_overlap);
		}
		if(PorQ == 2 && n_terms>1){
			double *v_overlap = init_array(bstringcount +1);
			get_overlapP(&(c0[0][0]),v[0],&(P[0][(n_terms-1)*astringcount]),&(Q[0][0]), &(P[0][0]), v_overlap, n_terms,alphae, betae, nmo);
			double final_overlap = C_DDOT(bstringcount+1,v[0],1,v_overlap,1);
			fprintf(outfile,"Final Eigenvector Overlap = %lf\n",C_DDOT(bstringcount+1,v[0],1,v_overlap,1));
			if(fabs(1.0-final_overlap)>10E-6){
				printf("Wavefunction is not normalized: norm = %lf\n",final_overlap);
			}
			free(v_overlap);
		}*/
			free_block(v);
				
	   }
	if(n_terms>1){
		free_block(eig_overlap);
	}
	free_block(b);
	free_block(sigma);
	free_block(G);
	free_block(alpha);
	free(lambda);
	free(old_lambda);
			return 0;

}

