/* This routine builds up a factored wavefunction one outer-product term at a time, optimizing the P and Q tables for term until the energies are self-consistent
 *
 * Using full length sigma builds
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

#include "lib/rfci_method.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;
int rfci_slow_method(mol_constant *mol, char * output, int print){
	//only solving the ground state for now
	int state = 0;
	//output file
	FILE * outfile = fopen(output,"a");
	//read settings file
	struct rfci_settings *rfci_set = read_rfci_settings();
	//initialize structs and open output files
	struct RFCI_wfn *wfn = allocate_memory(rfci_set,mol);
	wfn->max_dav_iters = *rfci_set->RFCI_MAX_DAV_ITERS;
	double rfci_cutoff = *rfci_set->RFCI_CUTOFF; 
	char Ptable_name[200];
	sprintf(Ptable_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/RFCI/tables/P_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
	FILE *Ptable_file = fopen(Ptable_name,"w");
	if(Ptable_file!=NULL){
		printf("Printing P tables to %s\n",Ptable_name);
	}
	else{
		printf("Could NOT print P tables to %s\n",Ptable_name);
	}

	char Qtable_name[200];
	sprintf(Qtable_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/RFCI/tables/Q_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
	FILE *Qtable_file = fopen(Qtable_name,"w");
	if(Qtable_file!=NULL){
		printf("Printing Q tables to %s\n",Qtable_name);
	}
	else{
		printf("Could NOT print Q tables to %s\n",Qtable_name);
	}

	char wfn_name[200];
	sprintf(wfn_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/RFCI/wfn/RFCIwfn_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
	FILE *wfn_file = fopen(wfn_name,"w");

	char energy_name[200];
	sprintf(energy_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/RFCI/energy/RFCIenergy_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
	FILE *energy_file = fopen(energy_name,"w");
	if(energy_file!=NULL){
		printf("Printing energies to %s\n",energy_name);
	}
	else{
		printf("Could NOT print energies tables to %s\n",energy_name);
	}
	//useful constants
	int astringcount = mol->astringcount; 
	int bstringcount = mol->bstringcount;
	int N = astringcount*bstringcount; 

	int wfn_terms, dav_iters = 0;
	int max_wfn_terms = *rfci_set->RFCI_MAX_WFN_ITERS;
	if(max_wfn_terms>astringcount*bstringcount){
		max_wfn_terms = astringcount*bstringcount;
	}
	int roots = *mol->roots;
	//verbose printing
	int debug_print = 1;

	int i,j,k,l;
	//difference between Davidson energies for a fixed number of wfn terms
	double diff = 1.0;
	//difference between Davidson energies for different numbers of wfn terms (finite difference derivative, kind of)
	double diff_of_diff  = 1.0;
	double last_energy;
	int num_wfn_terms;
	printf("max wfn terms = %d\n",max_wfn_terms);
	printf("\n     ======RFCI SLOW METHOD=====\n");
	wfn->c0[0][0] = 1.0;
	int guess_switch = 0;
	int rebuildQ = 1;

	//TRY TO SVD THE FCI ANSWER FOR BETTER INITIAL GUESS (ONLY RUN IF FCI HAS BEEN RUN ALREADY!)
	int svd_fci = 0;
	double *C_fci;
	double **U; 
	double **VT;
	if(svd_fci){
		char C_fci_name[200];
		sprintf(C_fci_name,"/home/settje/RFCI/molecule/%s/basis/%s/data/FCI/wfn/FCIwfn_%s_%s_%lf",mol->molname,mol->basisname,mol->molname,mol->basisname,*mol->rxn_coord);
		ifstream C_fci_file;
		C_fci_file.open(C_fci_name);
		if(C_fci_file.is_open()){ //read FCI wfn from file
			C_fci = init_array(astringcount*bstringcount*astringcount*bstringcount);
			printf("Reading FCI wavefunction from file: %s\n",C_fci_name);
			double buff;
			int C_ind = 0;
			while(C_fci_file >> buff){
				C_fci[C_ind] = buff;
				C_ind++;
			}
		}
		double *S = init_array(astringcount);
		U = block_matrix(astringcount,bstringcount);
		VT = block_matrix(bstringcount,astringcount);
		double *work = init_array(100*astringcount);
		int *iwork = init_int_array(100*astringcount);
		int info = 0;
		//SVD
		C_DGESDD('A',astringcount,bstringcount,C_fci,bstringcount,S,U[0],bstringcount,VT[0],astringcount,work,100*astringcount,iwork,info);
		fprintf(outfile,"U\n");
		print_mat(U,astringcount,bstringcount,outfile);
		fprintf(outfile,"V^T\n");
		print_mat(VT,astringcount,bstringcount,outfile);
		for(i=0;i<astringcount;i++){
			printf("S(%d) = %lf\n",i,S[i]);
		}
		for(i=0;i<astringcount;i++){
			wfn->P[0][i] = U[0][i];
		}
		normalize(wfn->P[0],astringcount);
		for(i=0;i<bstringcount;i++){
			wfn->Q[0][i] = VT[0][i];
		}
		normalize(wfn->Q[0],bstringcount);
	}	
	//MAJOR LOOP
	for(num_wfn_terms=0;num_wfn_terms<max_wfn_terms;num_wfn_terms++){
	//initialize guesses for current iteration
		wfn->c0[0][num_wfn_terms+1] = 1.0;
		if(num_wfn_terms==0 && guess_switch == 0 && svd_fci == 0){
			for(i=0;i<roots;i++){
				for(j=0;j<astringcount;j++){
					wfn->P[i][num_wfn_terms*astringcount+j] = 1/sqrt(astringcount);
				}
				for(j=0;j<num_wfn_terms;j++){
					double proj = C_DDOT(astringcount,&(wfn->P[i][num_wfn_terms*astringcount]),1,&(wfn->P[i][j*astringcount]),1);
					C_DAXPY(astringcount,-proj,&(wfn->P[i][j*astringcount]),1,&(wfn->P[i][num_wfn_terms*astringcount]),1);
				}
				normalize(&(wfn->P[i][num_wfn_terms*astringcount]),astringcount);
			}
		}
		if(rebuildQ && svd_fci == 0){
			for(i=0;i<roots;i++){
				for(j=0;j<bstringcount;j++){
					if(svd_fci){wfn->Q[i][num_wfn_terms*bstringcount+j] = VT[i][j];}
					else{
						wfn->Q[i][num_wfn_terms*bstringcount+j] = 1/sqrt(bstringcount);
					}
				}
				for(j=0;j<num_wfn_terms;j++){
					double proj = C_DDOT(bstringcount,&(wfn->Q[i][num_wfn_terms*bstringcount]),1,&(wfn->Q[i][j*bstringcount]),1);
					C_DAXPY(bstringcount,-proj,&(wfn->Q[i][j*bstringcount]),1,&(wfn->Q[i][num_wfn_terms*bstringcount]),1);
				}
				normalize(&(wfn->Q[i][num_wfn_terms*bstringcount]),bstringcount);
			}
		}
		else{
			rebuildQ = 1;
			svd_fci = 0;
		}
		
		printf("\n           NO. WFN TERMS = %d\n",num_wfn_terms+1);
		
		

		printf("It. Root           Ep                Eq            |Ep-Eq|\n");
		printf("-----------------------------------------------------------------------\n");
		dav_iters = 0;
		diff = 1.0;
		while( diff > rfci_cutoff && dav_iters <100){
			rfci_slow_davidson(1,wfn,mol,output,print);
			rfci_slow_davidson(2,wfn,mol,output,print);
			diff = wfn->Ep[0]-wfn->Eq[0]; //difference between energies calculated this iteration
			for(i=0;i<roots;i++){
				printf("%2d %2d %22.16lf %19.16lf %19.16lf\n",dav_iters,i,wfn->Ep[0],wfn->Eq[0],diff);
			}

			//FOR DEBUGGING
			if(num_wfn_terms>1 ){
				//return 0;
			}
			//END DEBUGGING
			dav_iters++;
		}
		if(energy_file!=NULL){
		for(i=0;i<roots;i++){
			printf("Printing energy %lf to %s\n",wfn->Eq[0],energy_name);
			fprintf(energy_file,"%20.16lf\n",wfn->Eq[0]);
		}
		}
		//PRINT TABLES AND ENERGIES
		if(Ptable_file!=NULL){
		for(i=0;i<astringcount;i++){
			fprintf(Ptable_file,"%20.16lf\n",wfn->P[0][num_wfn_terms*astringcount+i]);
		}
		}
		if(Qtable_file!=NULL){
		for(i=0;i<bstringcount;i++){
			fprintf(Qtable_file,"%20.16lf\n",wfn->Q[0][num_wfn_terms*astringcount+i]);
		}
		}

		//rescale wfn by optimum c0
		if(num_wfn_terms>0){
			printf("%lf ",wfn->c0[0][num_wfn_terms-1]);
			for(i=0;i<bstringcount;i++){
				printf("%lf ",wfn->Q[0][num_wfn_terms*bstringcount+i]);
			}
			printf("\n");
			for(i=0;i<roots;i++){
				C_DSCAL(num_wfn_terms*bstringcount,wfn->c0[i][num_wfn_terms-1],&wfn->Q[i][0],1);
				wfn->c0[i][num_wfn_terms-1] = 1.0;
			}
			
		}
		if(num_wfn_terms > 0 ){
			diff_of_diff = last_energy - wfn->Eq[0];
		}
		
		last_energy = wfn->Eq[0];
		wfn->n_terms++; //increment wfn terms in wfn struct too
	}

	if(num_wfn_terms<N){
		for(i=num_wfn_terms+1;i<N;i++){
			fprintf(energy_file,"%20.16lf\n",wfn->Eq[0]);
		}
	}
	free(C_fci);
	free_block(U);
	free_block(VT);
	if(Ptable_file!=NULL){
	fclose(Ptable_file);
	}
	if(Qtable_file!=NULL){
	fclose(Qtable_file);
	}
	if(wfn_file!=NULL){
	fclose(wfn_file);
	}
	if(energy_file!=NULL){
	fclose(energy_file);
	}
	
	fclose(outfile);
	free_rfci_settings(rfci_set);
	free_rfci_wfn(wfn);
	return 0;
}

