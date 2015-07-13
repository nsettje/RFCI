/* This routine builds up a factored wavefunction one outer-product term at a time, optimizing the P and Q tables for term until the energies are self-consistent
 *
 * Using factored length sigma builds
 *
 * output:
 * file containing optimum energy for wfn at each term
 * file containing final P tables
 * file containing final Q tables
 *
 * input:
 * mol struct of constants
 * output file name
 * print (switch for verbose terminal)
 */


#include "lib/RFCI.h"
#include "lib/rfci_wfn_memory.h"
#include "lib/rfci_settings.h"
//runs variational matrix decomposition algorithm by adding one alpha and one beta determinant at a time
int rfci_method(mol_constant *mol, char * output, int print){
	int state = 0;
	FILE * outfile = fopen(output,"w");
	struct rfci_settings *rfci_set = read_rfci_settings();
	struct RFCI_wfn *wfn = allocate_memory(rfci_set,mol);
	int i,j,k,l;
	
	//read in constants from the mol struct
	int nmo = *mol->nmo;
	int alphae = *mol->alphae;
	int betae = *mol->betae;

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	
	int wfn_terms, dav_iters;
	int max_wfn_terms = *rfci_set->RFCI_MAX_WFN_ITERS;
	int max_dav_iters = *rfci_set->RFCI_MAX_DAV_ITERS;

	//arrays to keep track of convergence in each Davidson iteration
	double *dav_Ep = init_array(max_dav_iters);
	double *dav_Eq = init_array(max_dav_iters);
	double *dav_diff = init_array(max_dav_iters);

	//number of Davidson roots
	int roots = *mol->roots;

	//converged tables
	double **Q = wfn->Q; 
	double **P = wfn->P; 

	//energy tolerance for convergence
	double Etol = *rfci_set->RFCI_CUTOFF;

	//best guesses for energy from P and Q iterations
	double *Ep = wfn->Ep; 
	double *Eq = wfn->Eq; 

	//normalization coefficients	
	double **c0 = wfn->c0; 

	double **curr_norm = block_matrix(max_dav_iters,4);



	//converged energies from all iterations	
	double *conv_energies = init_array(max_wfn_terms);

	//begin iterations, adding one new term as each converges
	for(wfn_terms=0;wfn_terms<max_wfn_terms;wfn_terms++){
		//all initializations for P and Q happen outside of the Davidson function
		//
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
			//rfci_davidson(1, roots, P, Q, wfn_terms+1, alphae, betae, nmo, mol->mo_OEI, mol->mo_OEIprime, mol->mo_TEI,c0,Ep);
			rfci_davidson(1,wfn,mol,output,1);
			
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
			//rfci_davidson(2, roots, P, Q, wfn_terms+1, alphae, betae, nmo, mol->mo_OEI, mol->mo_OEIprime, mol->mo_TEI,c0,Eq);
			rfci_davidson(2,wfn,mol,output,1);
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
		//FILE *energyout = fopen(energyoutfile,"a");
		//fprintf(energyout,"%20.16lf\n",conv_energies[wfn_terms]);
		//fclose(energyout);
		
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
	free_rfci_settings(rfci_set);
	free_rfci_wfn(wfn);
	return max_wfn_terms;
	free(conv_energies);

}

