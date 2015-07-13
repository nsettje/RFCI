#include "MOint.h"
#include "permute.h"
#include "slaterd.h"

/*Matrix Decomposition Full Configuration Interaction
 Based on Koch, Henrik, and Esper Dalgaard. "A variational matrix decomposition applied to full configuration-interaction calculations." Chemical physics letters 198.1 (1992): 51-58.

Most of the functionality is an implementation of equation 24 and its counterpart in beta
*/

double normalize(double *vec, int length);
double dot_vec(double *vec1, double *vec2, int length);
void add_scale_vec(double *vec1,double *vec2, double scale, int length);
double randouble();
void expand_factored_wfn(double *C, double c0, double *P, double *Q, int alphae, int betae, int nmo);
void project_expanded_sigmaP(double *sigma,double *proj_sigma, double *Q, int alphae, int betae, int nmo);
void project_expanded_sigmaQ(double *sigma,double *proj_sigma, double *Q, int alphae, int betae, int nmo);
void test_factored_sigma(int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI);
int variational_matrix_decomposition(int state, int nmo, int alphae, int betae, double nuc_rep_energy, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,char *energyoutfile, int print);
int davidson(int PorQ, int M, double **P, double **Q, int n_terms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI,double **c0, double *energy);
int davidP(int state, int M, double *total_energy,double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,double *c0,int use_guess,int print);
int davidQ(int state, int M,  double *total_energy,double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,double *c0,int use_guess,int print);
double get_energyP(int state, double **P, double **sigmaP, int n_Aterms, int alphae, int nmo);
double get_energyQ(int state, double **Q, double **sigmaQ, int n_Bterms, int betae, int nmo);
int get_factored_sigmaP(double scale, double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print);
void get_factored_sigmaQ(double scale, double *P,double *Q,double *Pprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print);
void get_factored_sigmaPAA(double *P,double *Q,double *Qprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaPBB(double *P,double *Q,double *Qprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaPAB(double *P,double *Q,double *Qprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQAA(double *P,double *Q,double *Pprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQBB(double *P,double *Q,double *Pprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQAB(double *P,double *Q,double *Pprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
int get_NO_augmented_sigmaP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *sigmaP, int n_terms,int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI);
int get_NO_augmented_sigmaQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *sigmaQ, int n_terms,int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI);
int get_overlapP_twoterm(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *overlapP, int n_terms,int alphae, int betae, int nmo);
double normalize_over_metricP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, int n_terms,int alphae, int betae, int nmo);
double normalize_over_metricQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, int n_terms,int alphae, int betae, int nmo);
int get_augmented_sigmaP(double *Pnew, double *Qnew, double *P, double *Q, double *sigmaP,  int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI);
int get_augmented_sigmaQ(double *Pnew, double *Qnew, double *P, double *Q, double *sigmaP,  int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI);
void test_augmented_sigma(int n_terms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI);
void test_NO_augmented_sigma(int n_terms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI);
extern "C" 
PsiReturnType variational_decomp_fci(Options& options)
{
    int print = options.get_int("PRINT");
/* INITIALIZATIONS */
	//pointers to important constants
	int *nmo_ = new int; //# molecular orbitals
	int *alphae_ = new int; //# alpha electrons
	int *betae_ = new int; //# beta electrons
	double *nuc_rep_energy_= new double;
	double *eSCF_ = new double;
	char *molname = new char[80];
	char *basisname = new char[80]; //molecule name and basis name

	//find important constants
	//function in MOint.h
	initialize_MO_constants(nmo_,alphae_,betae_,nuc_rep_energy_,eSCF_,molname,basisname);
	//recast pointers as ints
	int nmo = *(nmo_);
	int alphae = *alphae_;
	int betae = *betae_;
	double nuc_rep_energy= *(nuc_rep_energy_);
	double eSCF = *(eSCF_);
	//eSCF-=nuc_rep_energy;
	delete [] eSCF_;
	delete [] nmo_;
	delete [] alphae_;
	delete [] betae_;

	//initialize one-electron and two-electron integrals
	double **mo_OEI=block_matrix(nmo,nmo); //Pointer to OEI in MO basis
	double *mo_TEI=init_array(nmo*nmo*nmo*nmo);

	//transform MOs, storing in **OEI and *TEI arrays
	//function in MOint.h
	MO_transform(mo_OEI,mo_TEI,nmo);
	printf("Transformed MO integrals\n");

	//number of determinants
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);                                                          
	//full size of one dimension of NxN Hamiltonian  
	int N = astringcount*bstringcount;     
	printf("alpha strings = %d\nbeta strings = %d\nnmo = %d\n",astringcount,bstringcount,nmo);

	//option to build and diagonalize the full Hamiltonian to test convergence code
	double etot = 0.0; //store ground state energy
	int buildH = 0;
        if(buildH){                       
		double **H = block_matrix(astringcount*bstringcount,astringcount*bstringcount);
		build_full_Hamiltonian(H,mo_OEI,mo_TEI,alphae,betae,nmo);
		double *lambda = init_array(astringcount*bstringcount);
		double *work= init_array(10*astringcount*bstringcount);
		C_DSYEV('V','U',astringcount*bstringcount,&(H[0][0]),astringcount*bstringcount,lambda,work,10*astringcount*bstringcount);
		for(int i=0;i<astringcount*bstringcount;i++){
			fprintf(outfile,"%lf\n",lambda[i]);
		}
		etot = lambda[0];
		free(lambda);
		free(work);
		free_block(H);
	}
	
	//copy OEI to new array to be updated with TEI terms
	double **mo_OEIprime = block_matrix(nmo,nmo);
	C_DCOPY(nmo*nmo,&(mo_OEI[0][0]),1,&(mo_OEIprime[0][0]),1);	

	//add TEI to OEI to build matrix of hkl'
	for(int k =0;k<nmo;k++){
		for(int l=0;l<nmo;l++){
			for(int j =0;j<nmo;j++){
				mo_OEIprime[k][l]-=0.5*mo_TEI[k*nmo*nmo*nmo+j*nmo*nmo+j*nmo+l];
			}
		}
	}
	print_mat(mo_OEIprime,nmo,nmo,outfile);

	//functions to test sigma builds and augmented sigma builds with the c0 normalization parameter	
	int test_sig = 0;
	if(test_sig){
		test_factored_sigma(alphae, betae, nmo, mo_OEI, mo_OEIprime, mo_TEI);
		test_augmented_sigma(3, alphae, betae, nmo, mo_OEI,mo_OEIprime, mo_TEI);
		test_NO_augmented_sigma(5, alphae, betae, nmo, mo_OEI, mo_OEIprime, mo_TEI);
	}
	shared_ptr<Molecule> mol = Process::environment.molecule();
	//for diatomics, compute the bond length
	//double bondlength = fabs(2*mol->fz(1)/1.889726);
	double bondlength = fabs(mol->fz(1)-mol->fz(0))/1.889726;
	//bondlength*=bondlength;
	printf("r = %lf\n",bondlength);
	char energy_file[200];
	sprintf(energy_file,"/home/settje/psi_plugins/variational_decomp_fci/systems/%s/basis/%s/data/%s_%s_%lf",molname,basisname,molname,basisname,bondlength);
	printf("Printing to file: %s\n",energy_file);
	FILE *energy_out = fopen(energy_file,"w");
	fprintf(energy_out,"%lf\n",bondlength);
	fprintf(energy_out,"%20.16lf\n",nuc_rep_energy);
	fprintf(energy_out,"%20.16lf\n",eSCF);
	fclose(energy_out);	
	//back-and-forth Davidson
	printf("Beginning decompostion iterations\n");
	variational_matrix_decomposition(0,nmo,alphae,betae,nuc_rep_energy,mo_OEI,mo_OEIprime,mo_TEI,energy_file,0);
	energy_out = fopen(energy_file,"a");
	fprintf(energy_out,"%20.16lf\n",etot);
	fclose(energy_out);	

        if(buildH){       
		printf("FCI Energy = %20.14lf\n",etot);
		printf(" HF Energy = %20.14lf\n",eSCF);
	}                
	free(mo_TEI);
	free_block(mo_OEIprime);
	
/* END INITIALIZATIONS */
return Success;
}

//return random double
double randouble(){
	double F = (double)rand() / RAND_MAX;
	return F;
}

void project_expanded_sigmaP(double *sigma,double *proj_sigma, double *Q, int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int i,j,k,l;
	for(i=0;i<astringcount;i++){
		for(j=0;j<bstringcount;j++){
			proj_sigma[i] += sigma[i*astringcount+j]*Q[j]; 
		}
	}	
	
}

void project_expanded_sigmaQ(double *sigma,double *proj_sigma, double *P, int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int i,j,k,l;
	for(i=0;i<bstringcount;i++){
		for(j=0;j<astringcount;j++){
			proj_sigma[i] += sigma[j*astringcount+i]*P[j]; 
		}
	}	
	
}

void expand_factored_wfn(double *C, double c0, double *P, double *Q, int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int i,j,k,l;
	for(i=0;i<astringcount;i++){
		for(j=0;j<bstringcount;j++){
			C[i*astringcount+j] += c0*P[i]*Q[j];
		}
	}	
}





//runs variational matrix decomposition algorithm by adding one alpha and one beta determinant at a time
int variational_matrix_decomposition(int state, int nmo, int alphae, int betae, double nuc_rep_energy, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,char *energyoutfile, int print){


	int i,j,k,l;

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	
	int wfn_terms, dav_iters;
	int max_wfn_terms = bstringcount;
	int max_dav_iters = 20;

	//arrays to keep track of convergence in each Davidson iteration
	double *dav_Ep = init_array(max_dav_iters);
	double *dav_Eq = init_array(max_dav_iters);
	double *dav_diff = init_array(max_dav_iters);

	//number of Davidson roots
	int roots = 1;

	//converged tables
	double **Q = block_matrix(roots,max_wfn_terms*bstringcount);
	double **P = block_matrix(roots,max_wfn_terms*astringcount);

	//energy tolerance for convergence
	double Etol = 10E-6;

	//best guesses for energy from P and Q iterations
	double *Ep = init_array(roots);
	double *Eq = init_array(roots);

	//normalization coefficients	
	double **c0 = block_matrix(roots,max_wfn_terms+1);

	double **curr_norm = block_matrix(max_dav_iters,4);



	//converged energies from all iterations	
	double *conv_energies = init_array(max_wfn_terms);


	//begin iterations, adding one new term as each converges
	for(wfn_terms=0;wfn_terms<max_wfn_terms;wfn_terms++){
		//all initializations for P and Q happen outside of the Davidson function
		//choose homogeneous guess 
		for(i=0;i<bstringcount;i++){
			Q[state][wfn_terms*bstringcount+i]=1/sqrt(bstringcount);
		}
		for(i=0;i<wfn_terms;i++){
			for(j=0;j<i;j++){
				double proj = C_DDOT(bstringcount,&(Q[state][j*bstringcount]),1,&(Q[state][wfn_terms*bstringcount]),1);
				C_DAXPY(bstringcount,-proj,&(Q[state][j*bstringcount]),1,&(Q[state][wfn_terms*bstringcount]),1);
			}
			normalize(&(Q[state][bstringcount*wfn_terms]),bstringcount);
		}
		if(wfn_terms==0){
			for(i=0;i<astringcount;i++){
				P[state][wfn_terms*astringcount+i]=1/sqrt(astringcount);
			}
		}

		fprintf(outfile,"---------------------------WFN TERMS = %d---------------------------\n",wfn_terms);
		//print current tables
		fprintf(outfile, "Current P Tables\n");
		print_mat(P,1,(wfn_terms+1)*astringcount,outfile);
		fprintf(outfile, "Current Q Tables\n");
		print_mat(Q,1,(wfn_terms+1)*bstringcount,outfile);
		//difference between P and Q energies
		double curr_energy_diff = 1.0;
		//initialize the normalization constants to 1.0
		c0[state][wfn_terms] = 1.0;
		//iterate over P and Q, holding each fixed in turn until their energies are self consistent
		dav_iters = 0;	
		while(curr_energy_diff>Etol && dav_iters < max_dav_iters){
			//Q is fixed
			printf("WFN TERMS = %d (P Davidson Iteration = %d)\n",wfn_terms+1,dav_iters);
			printf("---------------------------------------------\n");
			davidson(1, roots, P, Q, wfn_terms+1, alphae, betae, nmo, mo_OEI, mo_OEIprime, mo_TEI,c0,Ep);
			
			//store ground state energy
			dav_Ep[dav_iters] = Ep[0];

			//print P table
			if(print>1){
				printf("P wavefunction\n");
				for(j=0;j<astringcount;j++){
					printf("%lf ",P[state][wfn_terms*astringcount+j]);
				}
				printf("\n");
			}
			printf("(P norm = %lf)\n",C_DDOT(astringcount,&(P[state][wfn_terms*astringcount]),1,&(P[state][wfn_terms*astringcount]),1));
							
			printf("\nWFN TERMS = %d (Q Davidson Iteration = %d)\n",wfn_terms+1,dav_iters);
			printf("---------------------------------------------\n");
		
			//P is fixed
			davidson(2, roots, P, Q, wfn_terms+1, alphae, betae, nmo, mo_OEI, mo_OEIprime, mo_TEI,c0,Eq);
			dav_Eq[dav_iters] = Eq[0];
			if(wfn_terms == 1 && dav_iters == 0){
				//break;
			}
			//compute difference between P and Q iterations
			curr_energy_diff = fabs(Ep[0]-Eq[0]);

			//store for later printing	
			dav_diff[dav_iters] = curr_energy_diff;

			//check whether energy is decreasing; print
			if(dav_iters>0 && dav_diff[dav_iters-1]<dav_diff[dav_iters]){
				for(i=0;i<dav_iters+1;i++){
					printf("%d %lf %lf %lf\n",i,dav_Ep[i],dav_Eq[i],dav_diff[i]);
				}
				printf("NOT VARIATIONAL!\n");
			}

			if(print>1){
				printf("Q wavefunction\n");
				for(j=0;j<bstringcount;j++){
					printf("%lf ",Q[state][wfn_terms*bstringcount+j]);
				}
				printf("\n");
			}

			//check convergence
			if(fabs(curr_energy_diff)<Etol){
				printf("Convergence Reached: Delta = %lf\n",curr_energy_diff);
				break;
			}
			dav_iters++;

		}//end of Davidson iterations

		//zero print array elements from last iteration
		for(i=dav_iters+1;i<max_dav_iters;i++){
			dav_Ep[i] = 0.0;
			dav_Eq[i] = 0.0;
			dav_diff[i] = 0.0;
		}

		//print convergence pattern for current iteration
		printf("Iter	Ep	Eq	deltaE	\n");
		for(i=0;i<dav_iters+1;i++){
			printf("%d %lf %lf %lf\n",i,dav_Ep[i],dav_Eq[i],dav_diff[i]);
		}
	
		//scale previous normalization coefficients	
		for(i=0;i<wfn_terms;i++){
			c0[0][i]*=c0[0][wfn_terms];
		}

		//normalize current P table and scale normalization coefficients by the normalization required for Q. 
		//This prevents overflow in the P and Q table elements.	
		if(wfn_terms > 0){
			c0[0][wfn_terms] = normalize(&(P[state][wfn_terms*astringcount]),astringcount);
			c0[0][wfn_terms] *= normalize(&(Q[state][wfn_terms*bstringcount]),bstringcount);
		}

		//print normalization coefficients
		printf("Normalization Coefficients\n");
		for(i=0;i<wfn_terms+1;i++){
			printf("%lf ",c0[0][i]);
		}
		printf("\n");

		//store energies from each matrix decomp iteration for later printing
		conv_energies[wfn_terms] = Eq[0];

		//print current energies to file
		FILE *energyout = fopen(energyoutfile,"a");
		fprintf(energyout,"%20.16lf\n",conv_energies[wfn_terms]);
		fclose(energyout);
		
		//check convergence
		if(dav_iters<max_dav_iters){
			printf("								|||| MATRIX DECOMP ITERATION %d CONVERGED IN %d DAVIDSON ITERATIONS ||||\n",wfn_terms,dav_iters);
		}
		else{
			printf("		++++ MATRIX DECOMP ITERATION %d DID NOT CONVERGE ++++\n",wfn_terms);
			printf("\nVARIATIONAL MATRIX DECOMPOSITION RESULTS\n");
			printf("Terms | Energy\n");
			printf("______________\n");
			for(i=0;i<max_wfn_terms;i++){
				printf("%5d | %lf\n",i,conv_energies[i]);
			}
			printf("\n");
			return 0;
		}

		//re-initialize convergence criterion
		curr_energy_diff = 1.0;
		
	}
	printf("\nVARIATIONAL MATRIX DECOMPOSITION RESULTS\n");
	printf("Terms | Energy\n");
	printf("______________\n");
	for(i=0;i<max_wfn_terms;i++){
		printf("%5d | %20.14lf\n",i,conv_energies[i]);
	}
	printf("\n");
	return max_wfn_terms;
	free(conv_energies);

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

//given a generalized eigenvalue problem AX = BXe, this function computes matrix elements of XBX in order to normalize vectors over the metric
double normalize_over_metricP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, int n_terms,int alphae, int betae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	double *overlap = init_array(astringcount+1);
	get_overlapP(c0,Pnew, Qnew, P, Q, overlap, n_terms,alphae, betae, nmo);
	double norm = C_DDOT(astringcount+1,Pnew,1,overlap,1);
	free(overlap);
	norm = 1/sqrt(norm);
	C_DSCAL(astringcount+1,norm,Pnew,1);
	return norm;
}
double normalize_over_metricQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, int n_terms,int alphae, int betae, int nmo){
	int bstringcount = nchoosek(nmo,betae);
	double *overlap = init_array(bstringcount+1);
//	get_overlapQ(c0,Pnew, Qnew, P, Q, overlap, n_terms,alphae, betae, nmo);
	get_overlapP(c0,Qnew, Pnew, Q, P, overlap, n_terms,alphae, betae, nmo);
	double norm = C_DDOT(bstringcount+1,Qnew,1,overlap,1);
	free(overlap);
	norm = 1/sqrt(norm);
	C_DSCAL(bstringcount+1,norm,Qnew,1);
	return norm;
}

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

int get_overlapP_twoterm(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *overlapP, int n_terms,int alphae, int betae, int nmo){
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	//build |psi_0> components for P and Q
	double *P0 = init_array(astringcount);
	double *Q0 = init_array(bstringcount);
	for(i=0;i<n_terms-1;i++){
		C_DAXPY(astringcount,c0[i],&P[i*astringcount],1,P0,1); //P0 = sum_i^n-1 c_i * P_i
		C_DAXPY(astringcount,1.0,&Q[i*bstringcount],1,Q0,1); //Q0 = sum_i^n-1 Q_i
	}
	//compute overlap[0] = c0 * <P0|P0> * <Q0|Q0> + <P0|P> * <Q0|Q>
	overlapP[0] = Pnew[0]*C_DDOT(astringcount,P0,1,P0,1)*C_DDOT(bstringcount,Q0,1,Q0,1); 
	overlapP[0] += C_DDOT(astringcount,P0,1,&Pnew[1],1)*C_DDOT(bstringcount,Q0,1,Qnew,1); 
	//compute remainder of overlap = c0 *<Q0|Q> * P0 + <Q|Q> * P 
	C_DAXPY(astringcount,Pnew[0]*C_DDOT(bstringcount,Q0,1,Qnew,1),P0,1,&overlapP[1],1);
	C_DAXPY(astringcount,C_DDOT(bstringcount,Qnew,1,Qnew,1),&Pnew[1],1,&overlapP[1],1);
	return 0;

}

int get_overlapP(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *overlapP, int n_terms,int alphae, int betae, int nmo){
	
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	for(i=0;i<astringcount;i++){
		overlapP[i] = 0.0;
	}
	for(i=0;i<n_terms-1;i++){
		for(j=0;j<n_terms-1;j++){
			overlapP[0]+=Pnew[0]*c0[i]*c0[j]*C_DDOT(astringcount,&(P[i*astringcount]),1,&(P[j*astringcount]),1)*C_DDOT(bstringcount,&(Q[i*bstringcount]),1,&(Q[j*bstringcount]),1); 
		}
		overlapP[0]+=c0[i]*C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[1],1)*C_DDOT(bstringcount,&(Q[i*bstringcount]),1,Qnew,1); 
		C_DAXPY(astringcount,c0[i]*C_DDOT(bstringcount,&(Q[i*astringcount]),1,Qnew,1),&(P[i*astringcount]),1,&overlapP[1],1);
	}
	C_DSCAL(astringcount,Pnew[0],&overlapP[1],1);
	C_DAXPY(astringcount,C_DDOT(bstringcount,Qnew,1,Qnew,1),&Pnew[1],1,&overlapP[1],1);
	return 0;

}
int get_overlapQ(double *c0,double *Pnew, double *Qnew, double *P, double *Q, double *overlapQ, int n_terms,int alphae, int betae, int nmo){
	
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	//build |psi_0> components for P and Q
	double *P0 = init_array(astringcount);
	double *Q0 = init_array(bstringcount);
	for(i=0;i<n_terms-1;i++){
		C_DAXPY(astringcount,c0[i],&P[i*astringcount],1,P0,1);
		C_DAXPY(astringcount,1.0,&Q[i*bstringcount],1,Q0,1);
	}
	//compute overlap[0]
	overlapQ[0] = Qnew[0]*C_DDOT(astringcount,P0,1,P0,1)*C_DDOT(bstringcount,Q0,1,Q0,1);
	overlapQ[0] += C_DDOT(astringcount,P0,1,Pnew,1)*C_DDOT(bstringcount,Q0,1,&Qnew[1],1);
	//compute remainder of overlap 
	C_DAXPY(bstringcount,Qnew[0]*C_DDOT(astringcount,P0,1,Pnew,1),Q0,1,&overlapQ[1],1);
	C_DAXPY(bstringcount,1.0,&Qnew[1],1,&overlapQ[1],1);
	return 1;

}

//compute augmented sigma^P (including first term for previous iteration's energy)
int get_augmented_sigmaP(double *Pnew, double *Qnew, double *P, double *Q, double *sigmaP, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI){
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	double norm_const = 0.0;
	for(i=0;i<n_Aterms-1;i++){
		norm_const += C_DDOT(bstringcount,&(Q[i*bstringcount]),1,&Qnew[0],1);
	}
	//norm_const-=2*norm_const*norm_const;
	norm_const*=norm_const;
	norm_const=1.0-norm_const;
	norm_const=1/norm_const;
	norm_const  = 1.0;
//	printf("P norm const  = %lf\n",sqrt(1/norm_const));	

	double E0 = 0.0;

	for(i=0;i<n_Aterms-1;i++){
		double qi = C_DDOT(bstringcount,&(Q[i*bstringcount]),1,&Qnew[0],1);
		get_factored_sigmaP(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Qnew[0],&sigmaP[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		double *sigPbuff = init_array(astringcount);
		get_factored_sigmaP(1.0,&(P[i*astringcount]),&Qnew[0],&Qnew[0],&sigPbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		C_DAXPY(astringcount,-qi,&sigPbuff[0],1,&sigmaP[1],1);
		free(sigPbuff);
		for(j=0;j<n_Aterms-1;j++){
			double qj = C_DDOT(bstringcount,&(Q[j*bstringcount]),1,&Qnew[0],1);
			sigPbuff = init_array(astringcount);
			get_factored_sigmaP(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&(Q[j*bstringcount]),&sigPbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			E0 += C_DDOT(astringcount,sigPbuff,1,&(P[j*astringcount]),1);
			free(sigPbuff);
			sigPbuff = init_array(astringcount);
			get_factored_sigmaP(1.0,&(P[i*astringcount]),&Qnew[0],&Qnew[0],&sigPbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			E0 += qi*qj*C_DDOT(astringcount,sigPbuff,1,&(P[j*astringcount]),1);
			free(sigPbuff);
			sigPbuff = init_array(astringcount);
			get_factored_sigmaP(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Qnew[0],&sigPbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			E0 -= (qi+qj)*C_DDOT(astringcount,&sigPbuff[0],1,&(P[j*astringcount]),1);
			free(sigPbuff);
		}
	}
	sigmaP[0] = norm_const*Pnew[0]*E0;
	sigmaP[0] += sqrt(norm_const)*C_DDOT(astringcount,&sigmaP[1],1,&Pnew[1],1);
	C_DSCAL(astringcount,sqrt(norm_const)*Pnew[0],&sigmaP[1],1);
	get_factored_sigmaP(1.0,&Pnew[1],&Qnew[0],&Qnew[0],&sigmaP[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	return 1;
}

//compute augmented sigma^Q (including first term for previous iteration's energy)
int get_augmented_sigmaQ(double *Pnew, double *Qnew, double *P, double *Q, double *sigmaQ, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI){
	
	int i,j,k,l;
	
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	double *ovp_test = init_array(astringcount);
	for(i=0;i<astringcount;i++){
		ovp_test[i] = P[i];	
	}
	double proj = C_DDOT(astringcount,ovp_test,1,&Pnew[0],1);
	C_DAXPY(astringcount,-proj,&Pnew[0],1,ovp_test,1);
	//printf("Q norm const exp = %lf\n",normalize(ovp_test,astringcount));
	free(ovp_test);
	double E0 = 0.0;

	double norm_const = 0.0;
	for(i=0;i<n_Bterms-1;i++){
		norm_const += C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[0],1);
	}
	norm_const*=norm_const;
	norm_const=1.0-norm_const;
	norm_const=1/norm_const;
	norm_const  = 1.0;
	//printf("Q norm const  = %lf\n",sqrt(1/norm_const));	
	for(i=0;i<n_Bterms-1;i++){
		double pi = C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[0],1);
		get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigmaQ[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		double *sigQbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&Pnew[0],&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		C_DAXPY(bstringcount,-pi,&sigQbuff[0],1,&sigmaQ[1],1);
		free(sigQbuff);
		for(j=0;j<n_Bterms-1;j++){
			double pj = C_DDOT(astringcount,&(P[j*astringcount]),1,&Pnew[0],1);
			sigQbuff = init_array(bstringcount);
			get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&(P[j*astringcount]),&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			E0 += C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[j*bstringcount]),1);
			free(sigQbuff);
			sigQbuff = init_array(bstringcount);
			get_factored_sigmaQ(1.0,&Pnew[0],&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			E0 += pi*pj*C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[j*bstringcount]),1);
			free(sigQbuff);
			sigQbuff = init_array(bstringcount);
			get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
			E0 -= (pi+pj)*C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[j*bstringcount]),1);
			free(sigQbuff);
		}
	}
	sigmaQ[0] = norm_const*Qnew[0]*E0;
	sigmaQ[0] += sqrt(norm_const)*C_DDOT(bstringcount,&sigmaQ[1],1,&Qnew[1],1);
	C_DSCAL(bstringcount,sqrt(norm_const)*Qnew[0],&sigmaQ[1],1);
	get_factored_sigmaQ(1.0,&Pnew[0],&Qnew[1],&Pnew[0],&sigmaQ[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	return 1;
	//<psi0|H|psi0>	
	/*for(i=0;i<n_Bterms-1;i++){
		for(j=0;j<n_Bterms-1;j++){
			double *sigQbuff = init_array(bstringcount);
			get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&(P[j*astringcount]),&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,3);
			double E0buff = C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[j*bstringcount]),1);
			printf("E0buff = %lf\n",E0buff);
			E0 += E0buff;
			printf("E0 = %lf\n",E0);
			free(sigQbuff);
		}
	}
	
	//<psi0|Pb><Pb|H|Pb'><Pb'|psi0>
	for(i=0;i<n_Bterms-1;i++){
		double P_overlap = C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[0],1);
		double *sigQbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		//-2<psi0|H|Pb><Pb|psi0>
		//<psi0|Pb><Pb|H|Pb'>
		E0 += (P_overlap*P_overlap-2*P_overlap)*C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[i*bstringcount]),1);
		C_DAXPY(bstringcount,-P_overlap,&sigQbuff[0],1,&sigmaQ[1],1);
		free(sigQbuff);
	}
	for(i=0;i<n_Bterms-1;i++){
		double P_overlap = C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[0],1);
		printf("P overlap = %lf\n",P_overlap);
		double *sigQbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		E0 += P_overlap*P_overlap*C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[i*bstringcount]),1);
			printf("E0 = %lf\n",E0);
		free(sigQbuff);
	}
	for(i=0;i<n_Bterms-1;i++){
		double P_overlap = C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[0],1);
		double *sigQbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		E0 += -2*P_overlap*C_DDOT(bstringcount,&sigQbuff[0],1,&(Q[i*bstringcount]),1);
			printf("E0 = %lf\n",E0);
		free(sigQbuff);
	}

	for(i=0;i<n_Bterms-1;i++){
		double P_overlap = C_DDOT(astringcount,&(P[i*astringcount]),1,&Pnew[0],1);
		double *sigQbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		C_DAXPY(bstringcount,-P_overlap,&sigQbuff[0],1,&sigmaQ[1],1);
		free(sigQbuff);
	}
	for(i=0;i<n_Bterms-1;i++){
		double *sigQbuff = init_array(bstringcount);
		get_factored_sigmaQ(1.0,&(P[i*astringcount]),&(Q[i*bstringcount]),&Pnew[0],&sigQbuff[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
		C_DAXPY(bstringcount,1.0,&sigQbuff[0],1,&sigmaQ[1],1);
		free(sigQbuff);
	}
	sigmaQ[0] = Qnew[0]*E0;
	sigmaQ[0] += C_DDOT(bstringcount,&sigmaQ[1],1,&Qnew[1],1);
	C_DSCAL(bstringcount,Qnew[0],&sigmaQ[1],1);
	get_factored_sigmaQAA(&Pnew[1],&Qnew[1],&Pnew[1],&sigmaQ[1],alphae,betae,nmo,mo_OEIprime,mo_TEI);
	return 1;*/
}

//compute sigma^P for a given state up to a given number of alpha and beta terms in the trial wavefunction, multiply each element of the final vector by the factor "scale"
int get_factored_sigmaP(double scale, double *P,double *Q,double *Qprime,double *sigmaP, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print){
	int astringcount = nchoosek(nmo,alphae);
	print = 3;
	get_factored_sigmaPAA(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>0){
		printf("sigmaP 1:     ");
		for(int i=0;i<astringcount;i++){
			//printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
	get_factored_sigmaPBB(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>1){
		printf("sigmaP 1+2:   ");
		for(int i=0;i<astringcount;i++){
			//printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}

	get_factored_sigmaPAB(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>2){
		printf("sigmaP 1+2+3: ");
		for(int i=0;i<astringcount;i++){
			//printf("%lf ",sigmaP[i]);
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


int davidson(int PorQ, int M, double **P, double **Q, int n_terms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI,double **c0, double *energy){
	int i, j, k, L, I;
	int min_pos, numf, iter, *conv, converged, maxdim, skip_check, built;
	int *small2big, init_dim;
	double *Adiag, **b,  **sigma, **G;
	double *lambda, **alpha, *f, *old_lambda,**eig_overlap;
	double norm, denom, diff;
	double BIGNUM = 10E100;
	int MAXIT = 1000;
	char vec_name;

	double cutoff = 10E-6;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N;
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
	maxdim = 20*M;
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
		sq_rsp(init_dim, init_dim, initG, lambda, 1, alpha, 1e-12);

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
					get_overlapP(&(c0[0][0]),&(b[i][0]), &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
					printf("Built sigma %d\n",built);
					built++;
				}
			}
			else{
				for(i=built;i<L;i++){
					get_NO_augmented_sigmaQ(&(c0[0][0]),&(P[0][(n_terms-1)*astringcount]),&(b[i][0]), &(P[0][0]), &(Q[0][0]), &(sigma[i][0]), n_terms,alphae, betae, nmo, mo_OEIprime, mo_TEI);
					get_overlapP(&(c0[0][0]), &(b[i][0]),&(P[0][(n_terms-1)*astringcount]),&(Q[0][0]),&(P[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
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
			C_DSYGVD(1,'V','U',L,&(G[0][0]),maxdim,&(G_overlap[0][0]),L,lambda,work,50*L,iwork,20*L);
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
			sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);
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
					get_overlapP(&(c0[0][0]),&(b[i][0]), &(Q[0][(n_terms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
					for(j=0;j<i;j++){
						double proj = C_DDOT(N,eig_overlap[i],1,b[j],1);
						C_DAXPY(N,-proj,b[j],1,b[i],1);
					}
					normalize_over_metricQ(&(c0[0][0]),&(Q[0][(n_terms-1)*bstringcount]), b[i],&(Q[0][0]), &(P[0][0]), n_terms,alphae, betae, nmo);
				}
			}
			else{
				for(i=1;i<new_dim;i++){	
					get_overlapP(&(c0[0][0]), &(b[i][0]),&(P[0][(n_terms-1)*astringcount]),&(Q[0][0]),&(P[0][0]), &(eig_overlap[i][0]), n_terms,alphae, betae, nmo);
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
			energy[i] = lambda[i];
			printf("%d %20.14lf\n",i,lambda[i]);
		}
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


//Davidson diagonalization for P table, holding Q fixed
int davidP(int state, int M, double *total_energy, double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI,double *c0,int use_guess,int print)
{
	//dummy variables
	int i, j, k, L, I;
	double minimum;	
	int min_pos, numf, iter, *conv, converged, maxdim, skip_check;
	int *small2big, init_dim;
	int smart_guess =1;
	double *Adiag, **b, **bnew, **sigma, **G;
	double *lambda, **alpha, **f, *lambda_old;
	double norm, denom, diff;
	double BIGNUM = 10E100;
	int add_flag;

	int guess_index;

	//maximum number of iterations
	int MAXIT = 100;
	if(n_Aterms>1){
		//MAXIT = 10;
	}
	
	//convergence criterion for the ENERGY
	double cutoff = 10E-9;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	//dimension includes c0 normalization parameter
	int N = astringcount+1;

	//maximum Davidson subspace dimension (M is the number of roots sought)
	if(n_Aterms == 1){
		maxdim = 25*M;
	}
	else{
		maxdim = (bstringcount-1) * M;
	}
	b = block_matrix(maxdim, N);  /* current set of guess vectors,
				   stored by row */
	bnew = block_matrix(2*M, N); /* guess vectors formed from old vectors,
				stored by row*/
	sigma = block_matrix(maxdim,N); /* sigma vectors, stored by row*/
	G = block_matrix(maxdim, maxdim); /* Davidson mini-Hamitonian */
	f = block_matrix(maxdim, N); /* residual eigenvectors, stored by row */
	alpha = block_matrix(maxdim, maxdim); /* eigenvectors of G */
	lambda = init_array(maxdim); /* eigenvalues of G */
	lambda_old = init_array(maxdim); /* approximate roots from previous
				      iteration */
	Adiag = init_array(N);
	

		if(n_Aterms>1){
		double norm_test = 0.0;
		norm_test+=c0[n_Aterms]*c0[n_Aterms]*C_DDOT(astringcount,&(P[state][(n_Aterms-1)*astringcount]),1,&(P[state][(n_Aterms-1)*astringcount]),1)*C_DDOT(bstringcount,&(Q[state][(n_Bterms-1)*bstringcount]),1,&(Q[state][(n_Bterms-1)*bstringcount]),1);
		norm_test+=c0[n_Aterms-1]*c0[n_Aterms-1];
		printf("c0 = %lf\n",c0[n_Aterms-1]);
		printf("total norm = %lf\n",norm_test);
		}
	if(smart_guess) { /* Use eigenvectors of a sub-matrix as initial guesses */
		//set initial dimension
		if(N > maxdim){
			init_dim = (maxdim-1)*M;
		}
		else{ 
			init_dim = 2*M;
		}
		init_dim = M;
		//if this switch is on, copy the previous best guess for current P
		//if(use_guess){
		//if(n_Aterms==1){
			C_DCOPY(astringcount,&(P[state][(n_Aterms-1)*astringcount]),1,&(b[0][1]),1);
	//	}
		C_DSCAL(astringcount,c0[n_Aterms],&(b[0][1]),1);
		//}
		//build sub-Hamiltonian element-by-element
		double **initG = block_matrix(init_dim,init_dim);
		for(i=0; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				initG[i][j] = build_single_Hamiltonian_element((n_Aterms-1)*astringcount+i,(n_Aterms-1)*astringcount+j,mo_OEI,mo_TEI,alphae,betae,nmo);
			}
		}
		//print_mat(initG,init_dim,init_dim,outfile);
		//C_DPBTRF('U',init_dim,init_dim-1,&(initG[0][0]),init_dim);
		fprintf(outfile,"init G P\n");
		print_mat(initG,init_dim,init_dim,outfile);

		//diagonalize sub-Hamiltonian
		sq_rsp(init_dim, init_dim, initG, lambda, 1, alpha, 1e-12);
		print_mat(alpha,init_dim,init_dim,outfile);
		
		//set the first element of b (c0) to 1, then normalize
		guess_index = 0;
		if(use_guess){
			guess_index = 1;
		}
		if(n_Aterms>1){
			b[0][0] = c0[n_Aterms-1];//sqrt(1.0-cP_scale[0][n_Aterms-1]*cP_scale[0][n_Aterms-1]);
			//normalize(b[0],N);
		}
		/*for(i=guess_index; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				b[i][j+1] = alpha[j][i-guess_index];
			}
			//b[i][0]=1.0;
			normalize(b[i],N);
		}*/
		
		int term_switch = 0;
		if(n_Aterms == 1){
			term_switch = 1;
		}
	
		if(n_Aterms==1){	
			for(i=0;i<astringcount;i++){
				b[0][i+1] = 1/sqrt(astringcount);
			}
		}
		/*for(i=0;i<n_Aterms-use_guess;i++){
			double proj = C_DDOT(astringcount,&(b[0][1]),1,&(P[state][i*astringcount]),1);
			C_DAXPY(astringcount,-proj,&(P[state][i*astringcount]),1,&(b[0][1]),1);
		}
		normalize(b[0],N);*/
	
		for(i=1;i<init_dim;i++){
			if(n_Aterms>1){
				for(j=0;j<N;j++){
					b[i][j] = 1/sqrt(N);
				}
			}
			for(j=0;j<i;j++){
				double proj = C_DDOT(N,&(b[j][0]),1,&(b[i][0]),1);
				C_DAXPY(N,-proj,&(b[j][0]),1,&(b[i][0]),1);
				//printf("Orthogonalized b(%d) wrt b(%d)\n",i,j);
			}
			/*for(j=0;j<n_Aterms-use_guess;j++){
				double proj = C_DDOT(astringcount,&(b[i][1]),1,&(P[state][j*astringcount]),1);
				C_DAXPY(astringcount,-proj,&(P[state][j*astringcount]),1,&(b[i][1]),1);
				//printf("Orthogonalized b(%d) wrt P(%d)\n",i,j);
			}*/
			normalize(b[i],N);
		}
		

		free_block(initG);
	}

	//set current dimension L to initial dimension
	L = init_dim;
	iter =0;
	converged = 0;
	conv = init_int_array(M); /* boolean array for convergence of each
			       root */
	//keep track of current number of sigma vectors in memory
	int built = 0;

	//build diagonal elements of H
	for(i=0;i<N;i++){
		Adiag[i] = build_single_Hamiltonian_element(i*astringcount+n_Bterms-1,i*astringcount+n_Bterms-1,mo_OEI,mo_TEI,alphae,betae,nmo);
	}


		printf("P wavefunction\n");
		for(i=0;i<astringcount+1;i++){
			printf("%lf ",b[0][i]);
		}
		printf("\n");
		printf("P wavefunction norm = %lf\n",C_DDOT(astringcount,b[0],1,b[0],1));

	//ITERATE
	double **energy_iters=block_matrix(MAXIT,2);
	int track_energy = 0;
	while(converged < M && iter < MAXIT) {
		//keep track of whether this iteration required subspace collapse
		skip_check = 0;
		if(print){
			 //printf("iter = %d\n", iter); 
		}
		//print guess vectors (b are guesses)
		fprintf(outfile,"B %d (L = %d)\n",iter,L);
		print_mat(b,L,N,outfile);
	
		//check orthonormality of b vectors
		int print_overlap = 1;
		if(print_overlap){
			for(i=0;i<n_Bterms;i++){
				for(j=0;j<n_Bterms;j++){
					//printf("P-P overlap (%d|%d) = %lf\n",i,j,C_DDOT(astringcount,&(P[state][i*astringcount]),1,&(P[state][j*astringcount]),1));
				}
			}
			for(i=0;i<L;i++){
				fprintf(outfile,"%lf\n",C_DDOT(astringcount,&(b[i][1]),1,&(P[state][(n_Aterms-1)*astringcount]),1));
				for(j=0;j<n_Aterms;j++){
					double bPov = C_DDOT(bstringcount,&(b[i][1]),1,&(P[state][j*astringcount]),1);
					if(fabs(bPov) > 10E-6 && i != 0){
						//printf("NONZERO OVERLAP b(%d) and P(%d) = %lf\n",i,j,bPov);
						//return 0;
					}
					//printf("b-P overlap (%d|%d) = %lf\n",i,j,C_DDOT(astringcount,&(b[i][1]),1,&(P[state][j*astringcount]),1));
				}
			}
			double **overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(b[0][0]), N, 0.0, &(overlap[0][0]), L);
			fprintf(outfile,"P overlap\n");
			print_mat(overlap,L,L,outfile);
			free_block(overlap);	
		}
		//build sigma vectors
		if(n_Aterms>1){
			for(i=built;i<L;i++){
				if(!get_augmented_sigmaP(&(b[i][0]), &(Q[state][(n_Bterms-1)*bstringcount]), &(P[0][0]), &(Q[0][0]), &(sigma[i][0]),  n_Aterms, n_Bterms, alphae, betae,  nmo, mo_OEIprime, mo_TEI)){return 0;}
				built++;
			}
		}
		else{
			for(i=built;i<L;i++){
				//printf("b %d: ",iter);
				for(j =0;j<astringcount;j++){
					//printf("%lf ",b[i][j+1]);
				}
				//printf("\n");
				if(!get_factored_sigmaP(1.0,&(b[i][1]),&(Q[state][0]),&(Q[state][0]),&(sigma[i][1]), alphae, betae, nmo, mo_OEIprime,  mo_TEI,0)){return 0;}
				built++;
			}
		}

		//print sigma vectors
		fprintf(outfile,"sigma P %d\n",iter);
		print_mat(sigma,built,N,outfile);
	
		/* form mini-matrix */
		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim);

		


		//print quadratic forms
		fprintf(outfile,"subHP %d\n",iter);
		print_mat(G,L,L,outfile);
	
		//check symmetry of mini-matrix
		double symm_check = 0.0;
		for(i=0;i<L;i++){
			for(j=0;j<L;j++){
				symm_check += fabs(G[i][j] - G[j][i]);
			}
		}
		if(fabs(symm_check)>10E-6){
			printf("P Davidson Sub-Hamiltonian is NOT SYMMETRIC!: symm = %lf\n",symm_check);
			return 0;
		}

		/* diagonalize mini-matrix */
		sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);
		fprintf(outfile,"P eigenvalues\n");
		for(j=0;j<M;j++){
			fprintf(outfile,"%lf\n",lambda[j]);
		}
		fprintf(outfile,"subHP eigenvectors\n");
		print_mat(alpha,L,L,outfile);

		/* form preconditioned residue vectors */
		for(k=0; k < M; k++) {//rows
			for(I=0; I < N; I++) { //cols
				f[k][I] = 0.0;
				for(i=0; i < L; i++) {
					f[k][I] += alpha[i][k] * (sigma[i][I] - lambda[k] * b[i][I]);
				}
				denom = lambda[k] - Adiag[I];
				if(fabs(denom) > 10e-6) {
					f[k][I] /= denom;
				}
				else{
					f[k][I] = 0.0;
				}
			}
		}

		/* normalize each residual */
		for(k=0; k < M; k++) {
			norm = 0.0;
			for(I=0; I < N; I++) {
				norm += f[k][I] * f[k][I];
			}
			norm = sqrt(norm);
			for(I=0; I < N; I++) {
				if(norm > 1e-6) {
					f[k][I] /= norm;
				}
				else {
					f[k][I] = 0.0;
				}
			}
		}

		/* schmidt orthogonalize the f[k] against the set of b[i] and add
		new vectors */
		add_flag = 0;
		for(k=0,numf=0; k < M; k++){
				for(I=0;I<L;I++){
					double proj = C_DDOT(N,&(b[I][0]),1,&(f[k][0]),1);
					C_DAXPY(N,-proj,&(b[I][0]),1,&(f[k][0]),1);
				}
				/*for(j=0;j<n_Aterms;j++){
					double proj = C_DDOT(astringcount,&(f[k][1]),1,&(P[state][j*astringcount]),1);
					C_DAXPY(astringcount,-proj,&(P[state][j*astringcount]),1,&(f[k][1]),1);
				}*/
				if(normalize(&(f[k][0]),N)>10E-6){
					for(i=0;i<N;i++){
						b[L][i] = f[k][i];
					}
					L++; numf++; 
				}
		}
		/* If L is close to maxdim, collapse to two guesses per root */
		if(maxdim - L < M) {
			if(print) {
				//printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
				//printf("Collapsing eigenvectors.\n");
				fprintf(outfile,"P COLLAPSE ==============\n");
			}
			for(i=0;i<L;i++){
				for(j=0;j<N;j++){
					sigma[i][j]=0.0;
				}
			}
			for(i=0; i < 2*M; i++) {
				memset((void *) bnew[i], 0, N*sizeof(double));
				for(j=0; j < L; j++) {
					for(k=0; k < N; k++) {
						bnew[i][k] += alpha[j][i] * b[j][k];
					}
				}
			}
			/* copy new vectors into place */
			for(i=0; i < 2*M; i++){ 
				for(k=0; k < N; k++){
					b[i][k] = bnew[i][k];
				}
			}
			if(n_Aterms==1){
				for(i=0;i<M;i++){
					b[i][0]=0.0;
					normalize(&(b[i][0]),N);
				}
			}
			for(i=0;i<2*M;i++){
				/*for(j=0;j<n_Aterms;j++){
					double proj = C_DDOT(astringcount,&(b[i][1]),1,&(P[state][j*astringcount]),1);
					C_DAXPY(astringcount,-proj,&(P[state][j*astringcount]),1,&(b[i][1]),1);
				}*/
			}
			normalize(&(b[0][0]),N);
			for(i=1;i<2*M;i++){
				for(j=0;j<i;j++){
					double proj = C_DDOT(N,&(b[i][0]),1,&(b[j][0]),1);
					C_DAXPY(N,-proj,&(b[j][0]),1,&(b[i][0]),1);
				}
				normalize(&(b[i][0]),N);
			}

					
			skip_check = 1;
			built = 0;
			L = 2*M;
		}

		/* check convergence on all roots */
		if(!skip_check) {
			energy_iters[track_energy][0] = lambda[0];
			if(track_energy>0){
				energy_iters[track_energy][1] = lambda[0] - energy_iters[track_energy-1][0];
			}
			converged = 0;
			zero_int_array(conv, M);
			if(print) {
				//printf("Root      Eigenvalue       Delta  Converged?\n");
				//printf("---- -------------------- ------- ----------\n");
			}
			for(k=0; k < M; k++) {
				diff = fabs(lambda[k] - lambda_old[k]);
				if(diff < cutoff) {
					conv[k] = 1;
					converged++;
				}
				lambda_old[k] = lambda[k];
				if(print) {
				//	printf("%3d  %20.14f %4.3e    %1s\n", k, lambda[k], diff,
					// conv[k] == 1 ? "Y" : "N");
				}
			}
			track_energy++;
		}

		//printf("c0 norm = %lf\n",C_DDOT(n_Aterms,cP_scale[0],1,cP_scale[0],1));
		iter++;
	}

	/* generate final eigenvalues and eigenvectors */
	if(converged == M) {
	//copy converged energies
	for(i=0;i<M;i++){
		total_energy[i] = lambda[i];
	}
	double **v = block_matrix(M,N);
		for(i=0; i < N; i++) {
			for(j=0; j < L; j++) {
				for(I=0; I < M; I++) {
					v[I][i] += alpha[j][I] * b[j][i];
				}
			}
		}
		//printf("eigenvector normalization = %20.16lf\n",normalize(&(v[0][0]),astringcount+1));
		for(i=0;i<M;i++){
			//curr_cP[i]  =  v[i][0];
		}
		if(print) printf("Davidson algorithm converged in %d iterations (including collapses).\n", iter);
		//copy newest converged vector
		for(I=0;I<astringcount;I++){
			P[state][(n_Aterms-1)*astringcount+I]=v[state][I+1];
		}
		double norm_const = 0.0;
		for(i=0;i<n_Aterms-1;i++){
			norm_const += C_DDOT(bstringcount,&(Q[state][i*bstringcount]),1,&Q[state][(n_Bterms-1)*bstringcount],1);
		}
		norm_const*=norm_const;
		norm_const=1.0-norm_const;
		norm_const=1/norm_const;
		printf("psi0' norm = %lf\n",1/norm_const);
		c0[n_Aterms-1] = v[0][0]*sqrt(norm_const);
		//c0[n_Aterms-1] = norm_const*v[0][0];
		
		//normalize newest converged vector
		c0[n_Aterms] = normalize(&(P[state][(n_Aterms-1)*astringcount]),astringcount);
		//printf("normalization factor = %lf\n",c0[n_Aterms]);
		printf("\n P eigenvector\n");
		for(i=0;i<astringcount+1;i++){
			printf("%lf ",v[0][i]);
		}
		printf("\n");
		printf("Q wavefunction\n");
		for(i=0;i<bstringcount;i++){
			printf("%lf ",Q[state][(n_Bterms-1)*bstringcount+i]);
		}
		printf("\n");

		fprintf(outfile,"||||||||||Converged eigenvectors P|||||||||||||\n");
	print_mat(v,M,N,outfile);
	fprintf(outfile,"|||||||||||||||||||||||||||||||||||||||||||||||\n");
	}
		printf("Iter		Ep		Delta\n");
		for(i=0;i<track_energy;i++){
			printf("%2d %20.14lf %20.14lf\n",i,energy_iters[i][0],energy_iters[i][1]);
		}	

  free(conv);
  free_block(b);
  free_block(bnew);
  free_block(sigma);
  free_block(G);
  free_block(f);
  free_block(alpha);
  free(lambda);
  free(lambda_old);

  return converged;
}


//same as davidP but it forms sigma^Q vectors
//see davidP for more extensive comments
int davidQ(int state, int M, double *total_energy, double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,double *c0,int use_guess,int print){
	int i, j, k, L, I;
	double minimum;
	int add_flag;
	int min_pos, numf, iter, *conv, converged, maxdim, skip_check;
	int *small2big, init_dim;
	int smart_guess =1;
	double *Adiag, **b, **bnew, **sigma, **G;
	double *lambda, **alpha, **f, *lambda_old;
	double norm, denom, diff;
	double BIGNUM = 10E100;
	int MAXIT = 100;
	double cutoff = 10E-9;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N = bstringcount+1;
	if(n_Aterms == 1){
		maxdim = 25*M;
	}
	else{
		maxdim = (astringcount-1)*M;
	}
	if(n_Bterms > 1){
		//MAXIT = 10;
	}
//	printf("normalization factor = %lf\n",c0[n_Bterms]);


	b = block_matrix(maxdim, N);  /* current set of guess vectors,
				   stored by row */
	bnew = block_matrix(2*M, N); /* guess vectors formed from old vectors,
				stored by row*/
	sigma = block_matrix(maxdim,N); /* sigma vectors, stored by row*/
	G = block_matrix(maxdim, maxdim); /* Davidson mini-Hamitonian */
	f = block_matrix(maxdim, N); /* residual eigenvectors, stored by row */
	alpha = block_matrix(maxdim, maxdim); /* eigenvectors of G */
	lambda = init_array(maxdim); /* eigenvalues of G */
	lambda_old = init_array(maxdim); /* approximate roots from previous
				      iteration */
		Adiag = init_array(N);
	if(n_Bterms == 1){
		smart_guess = 1;
		printf("Q Sub-Hamiltonian Guess\n");
	}

		if(n_Aterms>1){
		double norm_test = 0.0;
		norm_test+=c0[n_Bterms]*c0[n_Bterms]*C_DDOT(astringcount,&(P[state][(n_Aterms-1)*astringcount]),1,&(P[state][(n_Aterms-1)*astringcount]),1)*C_DDOT(bstringcount,&(Q[state][(n_Bterms-1)*bstringcount]),1,&(Q[state][(n_Bterms-1)*bstringcount]),1);
		norm_test+=c0[n_Aterms-1]*c0[n_Aterms-1];
		printf("c0 = %lf\n",c0[n_Aterms-1]);
		printf("total norm = %lf\n",norm_test);
		}

	
	if(smart_guess) { /* Use eigenvectors of a sub-matrix as initial guesses */
		if(N > maxdim){
			init_dim = (maxdim-1)*M;
		}
		else{ 
			init_dim = 2*M;
		}
		init_dim = M;
		if(use_guess){
			C_DCOPY(bstringcount,&(Q[state][(n_Bterms-1)*bstringcount]),1,&(b[0][1]),1);	
			C_DSCAL(astringcount,/*sqrt(1.0-c0[n_Aterms-1]*c0[n_Aterms-1])*/c0[n_Bterms],&(b[0][1]),1);
		}
		//build sub-Hamiltonian element-by-element
		double **initG = block_matrix(init_dim,init_dim);
			for(i=0; i < init_dim; i++) {
				for(j=0; j < init_dim; j++){
					initG[i][j] = build_single_Hamiltonian_element(i*astringcount+n_Bterms-1,j*astringcount+n_Bterms-1,mo_OEI,mo_TEI,alphae,betae,nmo);
				}
			}
		
		fprintf(outfile,"init G Q\n");
		print_mat(initG,init_dim,init_dim,outfile);
		
		//diagonalize sub-Hamiltonian
		sq_rsp(init_dim, init_dim, initG, lambda, 1, alpha, 1e-12);
		print_mat(alpha,init_dim,init_dim,outfile);
		
		//set c0 = 1 and normalize b vectors
		int guess_index=0;
		if(use_guess){
			guess_index = 1;
			//printf("Using previous Q table guess\n");
		}
		if(n_Aterms>1){
			b[0][0] = c0[n_Bterms-1];
			//normalize(b[0],N);
		}
		for(i=guess_index; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				//b[i][j+1] = alpha[j][i-guess_index];
			}
			//b[i][0]=1.0;
			//normalize(b[i],N);
		}
		print_mat(b,init_dim,N,outfile);
		printf("Q wavefunction\n");
		for(i=0;i<bstringcount+1;i++){
			printf("%lf ",b[0][i]);
		}
		printf("\n");
	
		int term_switch = 1;
		if(n_Bterms == 1){
			term_switch = 0;
		}

		/*for(i=0;i<n_Bterms-use_guess;i++){
			double proj = C_DDOT(bstringcount,&(b[0][1]),1,&(Q[state][i*bstringcount]),1);
			C_DAXPY(bstringcount,-proj,&(Q[state][i*bstringcount]),1,&(b[0][1]),1);
		}
		normalize(b[0],N);*/

		for(i=1;i<init_dim;i++){
			if(n_Bterms>1){
				for(j=0;j<N;j++){
					b[i][j] = 1/sqrt(N);
				}
			}
			for(j=0;j<i;j++){
				double proj = C_DDOT(N,&(b[j][0]),1,&(b[i][0]),1);
				C_DAXPY(N,-proj,&(b[j][0]),1,&(b[i][0]),1);
			}
			/*for(j=0;j<n_Bterms-use_guess;j++){
				double proj = C_DDOT(bstringcount,&(b[i][1]),1,&(Q[state][j*bstringcount]),1);
				C_DAXPY(bstringcount,-proj,&(Q[state][j*bstringcount]),1,&(b[i][1]),1);
			}*/
			normalize(b[i],N);
		}
			double **overlap = block_matrix(init_dim,init_dim);
			C_DGEMM('n','t', init_dim, init_dim, N, 1.0, &(b[0][0]), N,&(b[0][0]), N, 0.0, &(overlap[0][0]), init_dim);
			fprintf(outfile,"Q overlap\n");
			print_mat(overlap,init_dim,init_dim,outfile);
			free_block(overlap);	
			print_mat(Q,1,n_Bterms*bstringcount,outfile);

			

		free_block(initG);
	}

	L = init_dim;
	iter =0;
	converged = 0;
	conv = init_int_array(M); /* boolean array for convergence of each
			       root */
	int built = 0;

	//diagonal elements of H
	for(i=0;i<bstringcount;i++){
		Adiag[i] = build_single_Hamiltonian_element(astringcount*(n_Aterms-1)+i,astringcount*(n_Aterms-1)+i,mo_OEI,mo_TEI,alphae,betae,nmo);
	}

	//ITERATE
	double **energy_iters=block_matrix(MAXIT,2);
	int track_energy = 0;
	while(converged < M && iter < MAXIT) {

		skip_check = 0;
		if(print){
			 //printf("\niter = %d\n", iter); 
		}
		fprintf(outfile,"B Q %d\n",iter);
		print_mat(b,L,N,outfile);
		//check orthonormality of b
		int print_overlap = 1;
		if(print_overlap){
			for(i=0;i<n_Bterms;i++){
				for(j=0;j<n_Bterms;j++){
				//	printf("Q-Q overlap (%d|%d) = %lf\n",i,j,C_DDOT(bstringcount,&(Q[state][i*bstringcount]),1,&(Q[state][j*bstringcount]),1));
				}
			}
			for(i=0;i<L;i++){
				fprintf(outfile,"%lf\n",C_DDOT(bstringcount,&(b[i][1]),1,&(Q[state][(n_Bterms-1)*bstringcount]),1));
				for(j=0;j<n_Bterms;j++){
					double bQov = C_DDOT(bstringcount,&(b[i][1]),1,&(Q[state][j*bstringcount]),1);
					if(fabs(bQov) > 10E-6 && i !=0){
					/*	printf("\nNONZERO OVERLAP b(%d) and Q(%d) = %lf\n",i,j,bQov);
						printf("b(%d): ",i);
						for(k=0;k<bstringcount;k++){
							printf("%lf ",b[i][k+1]);
						}
						printf("\n");
						printf("Q(%d): ",j);
						for(k=0;k<bstringcount;k++){
							printf("%lf ",Q[state][j*bstringcount+k]);
						}
						printf("\n");*/
						//return 0;
					}
				}
			}
			double **overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(b[0][0]), N, 0.0, &(overlap[0][0]), L);
			fprintf(outfile,"Q overlap\n");
			print_mat(overlap,L,L,outfile);
			free_block(overlap);	
		}
		//build sigma vectors
		if(n_Bterms>1){
			for(i=built;i<L;i++){
				get_augmented_sigmaQ(&(P[state][(n_Aterms-1)*astringcount]), &(b[i][0]), &(P[0][0]), &(Q[0][0]), &(sigma[i][0]),  n_Aterms, n_Bterms, alphae, betae,  nmo, mo_OEIprime, mo_TEI);
				//get_augmented_sigmaP(&(b[i][0]),&(P[state][(n_Aterms-1)*astringcount]), &(Q[0][0]), &(P[0][0]), &(sigma[i][0]),  n_Aterms, n_Bterms, alphae, betae,  nmo, mo_OEIprime, mo_TEI);
				built++;
			}
		}
		else{
			for(i=built;i<L;i++){
				//printf("b %d: ",iter);
				for(int kk =1;kk<bstringcount+1;kk++){
				//	printf("%lf ",b[i][kk]);
				}
				//printf("\n");
				double *sigref = &(sigma[i][1]);
				get_factored_sigmaQ(1.0,&(P[state][0]),&(b[i][1]),&(P[state][0]),sigref,  alphae, betae, nmo, mo_OEIprime,  mo_TEI,0);
				built++;
			}
		}
		fprintf(outfile,"SIGMA Q %d\n",iter);
		print_mat(sigma,built,N,outfile);
		/* form mini-matrix */
		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim);
		fprintf(outfile,"subHQ %d\n",iter);
		print_mat(G,L,L,outfile);
		//check symmetry of mini-matrix
		double symm_check = 0.0;
		for(i=0;i<L;i++){
			for(j=0;j<L;j++){
				symm_check += fabs(G[i][j] - G[j][i]);
			}
		}
		if(fabs(symm_check)>10E-6){
			printf("Q Davidson Sub-Hamiltonian is NOT SYMMETRIC!: symm = %lf\n",symm_check);
			return 0;
		}
		/* diagonalize mini-matrix */
		sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);
		fprintf(outfile,"Q eigenvalues\n");
		for(int kk=0;kk<M;kk++){
			fprintf(outfile,"%lf\n",lambda[kk]);
		}
		fprintf(outfile,"subHQ eigenvectors\n");
		print_mat(alpha,L,L,outfile);
		/* form preconditioned residue vectors */
		for(k=0; k < M; k++) {//rows
			for(I=0; I < N; I++) { //cols
				f[k][I] = 0.0;
				for(i=0; i < L; i++) {
					f[k][I] += alpha[i][k] * (sigma[i][I] - lambda[k] * b[i][I]);
				}
				denom = lambda[k] - Adiag[I];
				if(fabs(denom) > 1e-6) {
					f[k][I] /= denom;
				}
				else{
					f[k][I] = 0.0;
				}
			}
		}

		/* normalize each residual */
		for(k=0; k < M; k++) {
			norm = 0.0;
			for(I=0; I < N; I++) {
				norm += f[k][I] * f[k][I];
			}
			norm = sqrt(norm);
			for(I=0; I < N; I++) {
				if(norm > 1e-6) {
					f[k][I] /= norm;
				}
				else {
					f[k][I] = 0.0;
				}
			}
		}

		/* schmidt orthogonalize the f[k] against the set of b[i] and add
		new vectors */
		for(k=0,numf=0; k < M; k++){
			//if(schmidt_add(b, L, N, f[k])) { 
				for(I=0;I<L;I++){
					double proj = C_DDOT(N,&(b[I][0]),1,&(f[k][0]),1);
					C_DAXPY(N,-proj,&(b[I][0]),1,&(f[k][0]),1);
				}
				/*for(j=0;j<n_Bterms;j++){
					double proj = C_DDOT(bstringcount,&(f[k][1]),1,&(Q[state][j*bstringcount]),1);
					C_DAXPY(bstringcount,-proj,&(Q[state][j*bstringcount]),1,&(f[k][1]),1);
				}*/
				if(normalize(&(f[k][0]),N)>10E-6){
				//normalize(&(f[k][0]),N);
					for(i=0;i<N;i++){
						b[L][i] = f[k][i];
					}
					L++; numf++; 
					add_flag = 1;
				}
				else{
					add_flag = 0;
				}
			//}
		}
		
		/* If L is close to maxdim, collapse to one guess per root */
		if(maxdim - L < M) {
			if(print) {
				//printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
				//printf("Collapsing eigenvectors.\n");
				fprintf(outfile,"Q COLLAPSE ==============\n");
			}
			for(i=0;i<L;i++){
				for(j=0;j<N;j++){
					sigma[i][j]=0.0;
				}
			}
			for(i=0; i < 2*M; i++) {
				memset((void *) bnew[i], 0, N*sizeof(double));
				for(j=0; j < L; j++) {
					for(k=0; k < N; k++) {
						bnew[i][k] += alpha[j][i] * b[j][k];
					}
				}
			}
			/* copy new vectors into place */
			for(i=0; i < 2*M; i++){ 
				for(k=0; k < N; k++){
					b[i][k] = bnew[i][k];
				}
			}
			if(n_Bterms==1){
				for(i=0;i<M;i++){
					b[i][0]=0.0;
					normalize(&(b[i][0]),N);
				}
			}
			for(i=0;i<2*M;i++){
				/*for(j=0;j<n_Bterms;j++){
					double proj = C_DDOT(bstringcount,&(b[i][1]),1,&(Q[state][j*bstringcount]),1);
					C_DAXPY(bstringcount,-proj,&(Q[state][j*bstringcount]),1,&(b[i][1]),1);
				}*/
			}
			normalize(&(b[0][0]),N);
			for(i=1;i<2*M;i++){
				for(j=0;j<i;j++){
					double proj = C_DDOT(N,&(b[i][0]),1,&(b[j][0]),1);
					C_DAXPY(N,-proj,&(b[j][0]),1,&(b[i][0]),1);
				}
				normalize(&(b[i][0]),N);
			}
			skip_check = 1;
			built = 0;
			L = 2*M;
		}

		/* check convergence on all roots */
		if(!skip_check) {
			energy_iters[track_energy][0] = lambda[0];
			if(track_energy>0){
				energy_iters[track_energy][1] = lambda[0] - energy_iters[track_energy-1][0];
			}
			converged = 0;
			zero_int_array(conv, M);
			if(print) {
				//printf("Root      Eigenvalue       Delta  Converged?\n");
				//printf("---- -------------------- ------- ----------\n");
			}
			for(k=0; k < M; k++) {
				diff = fabs(lambda[k] - lambda_old[k]);
				if(diff < cutoff) {
					conv[k] = 1;
					converged++;
				}
				
				lambda_old[k] = lambda[k];
				if(print) {
					//printf("%3d  %20.14f %4.3e    %1s\n", k, lambda[k], diff,
					 //conv[k] == 1 ? "Y" : "N");
				}
			}
			track_energy++;
		}

		iter++;
	}
	/* generate final eigenvalues and eigenvectors */
	if(converged == M ) {
	//copy converged energies
	for(i=0;i<M;i++){
		total_energy[i] = lambda[i];
	}
	double **v = block_matrix(N,M);
		for(i=0; i < M; i++) {
			for(j=0; j < L; j++) {
				for(I=0; I < N; I++) {
					v[I][i] += alpha[j][i] * b[j][I];
				}
			}
		}
		normalize(v[0],N);
		if(print) printf("Davidson algorithm converged in %d iterations (including collapses).\n", iter);
		for(I=0;I<bstringcount;I++){
			Q[state][(n_Bterms-1)*bstringcount+I]=v[I+1][state];
		}
		double norm_const = 0.0;
		for(i=0;i<n_Aterms-1;i++){
			norm_const += C_DDOT(astringcount,&(P[state][i*astringcount]),1,&P[state][(n_Aterms-1)*astringcount],1);
		}
		norm_const*=norm_const;
		norm_const=1.0-norm_const;
		norm_const=1/norm_const;
		printf("psi0' norm = %lf\n",1/norm_const);
		//c0[n_Aterms-1] = v[0][0];
		c0[n_Aterms-1] = v[0][0]/norm_const;
		//c0[n_Bterms-1] = v[0][0];
		c0[n_Bterms] = normalize(&(Q[state][(n_Bterms-1)*bstringcount]),bstringcount);
		printf("\nQ eigenvector\n");
		for(i=0;i<astringcount+1;i++){
			printf("%lf ",v[0][i]);
		}
		printf("\n");
		printf("P wavefunction\n");
		for(i=0;i<astringcount;i++){
			printf("%lf ",P[state][(n_Aterms-1)*astringcount+i]);
		}
		printf("\n");
		printf("P wavefunction norm = %lf\n",C_DDOT(astringcount,&P[state][(n_Aterms-1)*astringcount],1,&P[state][(n_Aterms-1)*astringcount],1));
		fprintf(outfile,"||||||||||||||||||||||Eigenvectors Q||||||||||||||||||||\n");
	print_mat(v,N,M,outfile);
	fprintf(outfile,"||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	}
		printf("Iter		Eq		Delta\n");
		for(i=0;i<track_energy;i++){
			printf("%2d %20.14lf %20.14lf\n",i,energy_iters[i][0],energy_iters[i][1]);
		}	

  free(conv);
  free_block(b);
  free_block(bnew);
  free_block(sigma);
  free_block(G);
  free_block(f);
  free_block(alpha);
  free(lambda);
  free(lambda_old);

  return converged;
}
void test_NO_augmented_sigma(int n_terms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI){

	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	double *P=init_array(n_terms*astringcount+1);
	double *Q=init_array(n_terms*bstringcount+1);
	//initialize P as random, with c0 stored in P[0]
	for(i=0;i<n_terms*astringcount+1;i++){
		P[i] = randouble();
	}
	//P[0] = 1.0;
	for(i=0;i<n_terms*bstringcount;i++){
		Q[i+1] = randouble();
	}
	for(i=1;i<n_terms;i++){
		normalize(&P[i*astringcount+1],astringcount);
		normalize(&Q[i*bstringcount+1],bstringcount);
	}
	normalize(&Q[1],bstringcount);	
	//normalize new P wrt the metric
	normalize(&P[0],astringcount+1);	
	printf("P: ");
	for(i=0;i<n_terms*astringcount+1;i++){
		printf("%lf ",P[i]);
	}
	printf("\n");
	printf("Q: ");
	for(i=0;i<n_terms*bstringcount+1;i++){
		printf("%lf ",Q[i]);
	}
	printf("\n");
	//populate the coefficients for each term as random
	double *c0 = init_array(n_terms-1);
	for(i=0;i<n_terms-1;i++){
		c0[i] = randouble();
	}
	//normalize the coefficients
	normalize(c0,n_terms-1);
	//compute sigmaP and sigmaQ
	double *sigmaP = init_array(astringcount+1);
	double *sigmaQ = init_array(bstringcount+1);
	get_NO_augmented_sigmaP(&(c0[0]),&(P[0]), &(Q[1]), &(P[astringcount+1]), &(Q[bstringcount+1]), &(sigmaP[0]), n_terms,alphae, betae, nmo, mo_OEIprime, mo_TEI);
	printf("sigmaP: ");
	for(i=0;i<astringcount+1;i++){
		printf("%lf ",sigmaP[i]);
	}
	printf("\n");
	get_NO_augmented_sigmaQ(&(c0[0]),&(Q[1]), &(P[0]), &(Q[bstringcount+1]), &(P[astringcount+1]), &(sigmaQ[0]), n_terms,alphae, betae, nmo, mo_OEIprime, mo_TEI);
	printf("sigmaQ: ");
	for(i=0;i<bstringcount+1;i++){
		printf("%lf ",sigmaQ[i]);
	}
	printf("\n");
	//build full Hamiltonian
	double **H = block_matrix(astringcount*bstringcount,astringcount*bstringcount);
	build_full_Hamiltonian(H,mo_OEI,mo_TEI,alphae,betae,nmo);
	
	//build "new" |P,Q>
	double *C = init_array(astringcount*bstringcount);
	expand_factored_wfn(C,1.0,&P[1],&Q[1],alphae,betae,nmo);
	
	//build sigma vector for "new" term
	double *sigma = init_array(astringcount*bstringcount);
	C_DGEMV('N',astringcount*bstringcount,astringcount*bstringcount,1.0,&(H[0][0]),astringcount*bstringcount,C,1,0,sigma,1);
	
	//build full C vector for each term of |psi_0>
	double *C0 = init_array(astringcount*bstringcount);
	
	//build sigma vector for |psi_0>
	double *sigma0 = init_array(astringcount*bstringcount);
	for(int i=1;i<n_terms;i++){
		expand_factored_wfn(C0,c0[i-1],&P[i*astringcount+1],&Q[i*bstringcount+1],alphae,betae,nmo);
	}
	C_DGEMV('N',astringcount*bstringcount,astringcount*bstringcount,1.0,&(H[0][0]),astringcount*bstringcount,C0,1,0,sigma0,1);
	double *exp_aug_sigmaP = init_array(astringcount+1);	
	exp_aug_sigmaP[0] = P[0]*C_DDOT(astringcount*bstringcount,C0,1,sigma0,1);
	exp_aug_sigmaP[0] += C_DDOT(astringcount*bstringcount,C0,1,sigma,1);
	project_expanded_sigmaP(sigma0,&exp_aug_sigmaP[1], &Q[1], alphae, betae, nmo);
	C_DSCAL(astringcount,P[0],&exp_aug_sigmaP[1],1);
	project_expanded_sigmaP(sigma,&exp_aug_sigmaP[1], &Q[1], alphae, betae, nmo);
	printf("expSIG: ");
	for(int i=0;i<astringcount+1;i++){
		printf("%lf ",exp_aug_sigmaP[i]);
	}
	printf("\n");

	//test overlaps
	double *P_overlap = init_array(astringcount+1);
	double *Q_overlap = init_array(bstringcount+1);
	get_overlapP(&(c0[0]),&(P[0]), &(Q[1]), &(P[astringcount+1]), &(Q[bstringcount+1]), &(P_overlap[0]), n_terms,alphae, betae, nmo);
	get_overlapQ(&(c0[0]),&(Q[1]), &(P[0]), &(Q[bstringcount+1]), &(P[astringcount+1]), &(Q_overlap[0]), n_terms,alphae, betae, nmo);
	printf("P overlap: ");
	for(i=0;i<astringcount+1;i++){
		printf("%lf ",P_overlap[i]);
	}
	printf("\n");
	printf("Q overlap: ");
	for(i=0;i<bstringcount+1;i++){
		printf("%lf ",Q_overlap[i]);
	}
	printf("\n");

	double *exp_overlap = init_array(astringcount+1);
	exp_overlap[0] = P[0]*C_DDOT(astringcount*bstringcount,C0,1,C0,1);
	exp_overlap[0] += C_DDOT(astringcount*bstringcount,C0,1,C,1);
	project_expanded_sigmaP(C0,&exp_overlap[1], &Q[1], alphae, betae, nmo);
	C_DSCAL(astringcount,P[0],&exp_overlap[1],1);
	project_expanded_sigmaP(C,&exp_overlap[1], &Q[1], alphae, betae, nmo);
	printf("  overlap: ");
	for(i=0;i<bstringcount+1;i++){
		printf("%lf ",exp_overlap[i]);
	}
	printf("\n");

	
}

//tests sigma builds augmented with c0 normalization parameter (eq. 18)
void test_augmented_sigma(int n_terms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI){
	srand(time(NULL));
	
	printf("--- test_augmented_sigma---\n");
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	double *P=init_array(n_terms*astringcount);
	double *sigmaP = init_array(astringcount+1);
	double *Q=init_array(n_terms*bstringcount);
	double *sigmaQ = init_array(bstringcount+1);
	double proj;
	//initialize P as random, with c0 stored in P[0]
	for(i=0;i<n_terms*astringcount;i++){
		P[i] = randouble();
	}
	for(i=0;i<n_terms;i++){
		normalize(&P[i*astringcount],astringcount);
	}
	
	//set P = Q
	for(int i=0;i<n_terms*astringcount;i++){
		Q[i] = P[i];
	}

	//populate vector of normalization coefficients and then normalize
	double *norm = init_array(n_terms);
	for(i=0;i<n_terms;i++){
		norm[i] = randouble();
	}
	normalize(norm,n_terms);
	for(i=0;i<n_terms;i++){
		C_DSCAL(astringcount,norm[i],&P[i*astringcount],1);
	}
	
	//normalize "old" wavefunction: |P|^2|Q|^2 = 1
	/*double old_norm = 0.0;
	for(int i=1;i<n_terms;i++){
		for(int j=1;j<n_terms;j++){
			old_norm += C_DDOT(astringcount,&P[i*astringcount],1,&P[j*astringcount],1)*C_DDOT(bstringcount,&Q[i*bstringcount],1,&Q[j*bstringcount],1);
		}
	}
	//enforce normalization through the P wavefunction
	for(int i=astringcount+1;i<(n_terms)*astringcount;i++){
		P[i]*=1/sqrt(old_norm);
	}*/
	//check normalization
	double old_norm = 0.0;
	for(int i=0;i<n_terms;i++){
		//for(int j=1;j<n_terms;j++){
			old_norm += C_DDOT(astringcount,&P[i*astringcount],1,&P[i*astringcount],1);
		//}
	}
	printf("old norm = %lf\n",old_norm);
	old_norm = 0.0;
	for(int i=1;i<n_terms;i++){
	//	for(int j=1;j<n_terms;j++){
			old_norm += C_DDOT(bstringcount,&Q[i*bstringcount],1,&Q[i*bstringcount],1);
	//	}
	}
	printf("old norm = %lf\n",old_norm);

	//normalize with respect to the metric: c0^2+|P|^2|Q|^2 = 1, where P and Q are from the "new" term	
	double norm_const = sqrt(1.0-P[0]*P[0]);
	for(int i=0;i<astringcount;i++){
		P[i+1]*=norm_const;
	}
	printf("P new norm = %lf\n",C_DDOT(astringcount+1,&P[0],1,&P[0],1));
	printf("P: ");
	for(int i=0;i<n_terms*astringcount+1;i++){
		printf("%lf ",P[i]);
	}
	printf("\n");
	printf("Q: ");
	for(int i=0;i<n_terms*bstringcount+1;i++){
		printf("%lf ",Q[i]);
	}
	printf("\n");
	
	
	//build augmented sigma^P
	get_augmented_sigmaP(&P[0], &Q[1], &P[astringcount+1], &Q[bstringcount+1], sigmaP, n_terms,n_terms, alphae, betae, nmo, mo_OEIprime, mo_TEI);
	printf("sigmaP: ");
	for(int i=0;i<astringcount+1;i++){
		printf("%lf ",sigmaP[i]);
	}
	printf("\n");
	
	//build augmented sigma^Q
	get_augmented_sigmaQ(&Q[1],&P[0],  &Q[bstringcount+1], &P[astringcount+1], sigmaQ, n_terms,n_terms, alphae, betae, nmo, mo_OEIprime, mo_TEI);
	printf("sigmaQ: ");
	for(int i=0;i<astringcount+1;i++){
		printf("%lf ",sigmaQ[i]);
	}
	printf("\n");
	
	//project out parts of "new" Q from "old" Q 
	for(int i=1;i<n_terms;i++){
		//newQ * oldQ
		proj = C_DDOT(bstringcount,&Q[0],1,&Q[i*bstringcount],1);
		//oldQ = oldQ - proj*newQ
		C_DAXPY(bstringcount,-proj,&Q[0],1,&Q[i*bstringcount],1);
	}
	for(int i=1;i<n_terms;i++){
		printf("Overlap Q(%d,%d) = %lf\n",0,i,C_DDOT(bstringcount,&Q[0],1,&Q[i*bstringcount],1));
	}

	//build full Hamiltonian
	double **H = block_matrix(astringcount*bstringcount,astringcount*bstringcount);
	build_full_Hamiltonian(H,mo_OEI,mo_TEI,alphae,betae,nmo);
	
	//build "new" |P,Q>
	double *C = init_array(astringcount*bstringcount);
	expand_factored_wfn(C,1.0,&P[1],&Q[1],alphae,betae,nmo);
	
	//build sigma vector for "new" term
	double *sigma = init_array(astringcount*bstringcount);
	C_DGEMV('N',astringcount*bstringcount,astringcount*bstringcount,1.0,&(H[0][0]),astringcount*bstringcount,C,1,0,sigma,1);
	
	//build full C vector for each term of |psi_0>
	double *C0 = init_array(astringcount*bstringcount);
	
	//build sigma vector for |psi_0>
	double *sigma0 = init_array(astringcount*bstringcount);
	for(int i=1;i<n_terms;i++){
		expand_factored_wfn(C0,1.0,&P[i*astringcount+1],&Q[i*bstringcount+1],alphae,betae,nmo);
	}
	C_DGEMV('N',astringcount*bstringcount,astringcount*bstringcount,1.0,&(H[0][0]),astringcount*bstringcount,C0,1,0,sigma0,1);
	double *exp_aug_sigmaP = init_array(astringcount+1);	
	exp_aug_sigmaP[0] = P[0]*C_DDOT(astringcount*bstringcount,C0,1,sigma0,1);
	exp_aug_sigmaP[0] += C_DDOT(astringcount*bstringcount,C0,1,sigma,1);
	project_expanded_sigmaP(sigma0,&exp_aug_sigmaP[1], &Q[1], alphae, betae, nmo);
	C_DSCAL(astringcount,P[0],&exp_aug_sigmaP[1],1);
	project_expanded_sigmaP(sigma,&exp_aug_sigmaP[1], &Q[1], alphae, betae, nmo);
	printf("expSIG: ");
	for(int i=0;i<astringcount+1;i++){
		printf("%lf ",exp_aug_sigmaP[i]);
	}
	printf("\n");
	/*double *lambda = init_array(astringcount*bstringcount);
	double *work= init_array(10*astringcount*bstringcount);
	C_DSYEV('V','U',astringcount*bstringcount,&(H[0][0]),astringcount*bstringcount,lambda,work,10*astringcount*bstringcount);
	for(int i=0;i<astringcount*bstringcount;i++){
		fprintf(outfile,"%lf\n",lambda[i]);
	}
	free(lambda);
	free(work);*/
	free_block(H);
}
//Computes <P',b|H|P,Q> and <a,Q'|H|P,Q> by multiplication with the explicit Hamiltonian and by get_factored_sigma and compares the results
//The purpose of using P, P', Q, and Q' is to ensure that the builds work with different bras and kets on either side of the matrix element
void test_factored_sigma(int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI){
	int i,j,k,l;
	printf("--- test_factored_sigma ---\n");
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	
	//build 2 random P and Q vectors, normalizing each separately
	double *P = init_array(2*astringcount);
	double *Q = init_array(2*bstringcount);
	for(i=0;i<2*astringcount;i++){
		P[i] = randouble();
	}
	normalize(&P[0],astringcount);
	normalize(&P[astringcount],astringcount);
	for(i=0;i<2*bstringcount;i++){
		Q[i] = randouble();
	}
	normalize(&Q[0],bstringcount);
	normalize(&Q[bstringcount],bstringcount);
	
	//build sigma P and sigma Q
	double *sigmaP = init_array(astringcount);
	double *sigmaQ = init_array(bstringcount);
	get_factored_sigmaP(1.0,&(P[0]),&(Q[0]),&Q[bstringcount],&sigmaP[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	get_factored_sigmaQ(1.0,&(P[0]),&(Q[0]),&P[astringcount],&sigmaQ[0],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	
	//expand P and Q into a full C vector
	double *C = init_array(astringcount*bstringcount);
	expand_factored_wfn(C,1.0,&P[0],&Q[0],alphae,betae,nmo);
	
	//build full Hamiltonian
	double **H = block_matrix(astringcount*bstringcount,astringcount*bstringcount);
	build_full_Hamiltonian(H,mo_OEI,mo_TEI,alphae,betae,nmo);

	//form full sigma vector
	double *sigma = init_array(astringcount*bstringcount);
	C_DGEMV('N',astringcount*bstringcount,astringcount*bstringcount,1.0,&(H[0][0]),astringcount*bstringcount,C,1,0,sigma,1);
	
	//project the full sigma onto P and Q
	double *proj_sigmaP = init_array(astringcount);
	double *proj_sigmaQ = init_array(bstringcount);
	project_expanded_sigmaP(sigma,proj_sigmaP, &Q[bstringcount], alphae, betae, nmo);
	project_expanded_sigmaQ(sigma,proj_sigmaQ, &P[astringcount], alphae, betae, nmo);

	//check agreement between the factored and exact sigma vectors
	C_DAXPY(astringcount,-1.0,proj_sigmaP,1,sigmaP,1);
	double P_error = C_DDOT(astringcount,sigmaP,1,sigmaP,1);
	if(fabs(P_error)<10E-6){
		printf("Factored SIGMA P and explicit SIGMA P agree.\n");
	}
	else{
		printf("Factored SIGMA P and explicit SIGMA P DO NOT agree.\n");
	}

	C_DAXPY(bstringcount,-1.0,proj_sigmaQ,1,sigmaQ,1);
	double Q_error = C_DDOT(bstringcount,sigmaQ,1,sigmaQ,1);
	if(fabs(Q_error)<10E-6){
		printf("Factored SIGMA Q and explicit SIGMA Q agree.\n");
	}
	else{
		printf("Factored SIGMA Q and explicit SIGMA Q DO NOT agree.\n");
	}

	free(P);
	free(Q);
	free(sigmaP);
	free(sigmaQ);
	free(C);
	free_block(H);
	free(sigma);
	free(proj_sigmaP);
	free(proj_sigmaQ);
}


}}



