/* Defines the struct of molecular constants mol_constant
 * Also defines functions for initializing the struct and reading constants from the output of PSI4
 *
 * Reads information from several files, so all initiations are paired with the files from which the values come
 *
 * output:
 * mol_constant struct with fields
 * 	molname
 * 	basisname
 * 	rxn_coord
 * 	method
 * 	roots
 * 	nmo
 * 	alphae
 * 	betae
 *	eSCF
 *	nuc_rep_energy
 *	mo_OEI (integrals)
 *	mo_TEI (integrals)
 *
 * input:
 * file name of PSI4 output with MO constants
 *
 */

#ifndef asdRFCI_io
#define asdRFCI_io
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "lib/memory.h"
#include "lib/mol_const.h"
#include "lib/permute.h"
using namespace std;

struct mol_constant *mem_input(char *input_name){
	//initialize
	struct mol_constant *input_params =  (mol_constant*) malloc(sizeof(struct mol_constant));
	input_params->molname = new char[100];
	input_params->basisname = new char[100];
	input_params->rxn_coord = new double;
	input_params->method = new char[100];
	input_params->roots = new int;
	FILE *input = fopen(input_name,"r");
	if(input == NULL){
		printf("File not found: %s\n",input_name);
	}
	//read
	if(input != NULL){
		char line[128];
		int inp_indx = 0;
		while(fgets(line,sizeof line,input)){
			if(line[0] != '#'){
				char *pos;
				if((pos=strchr(line,'\n')) != NULL){ *pos = '\0';}
				sprintf(line,"%s",line); //strips the '\n' off the end of the line
				if(inp_indx==0){
					sprintf(input_params->molname,"%s",line);
				}
				if(inp_indx==1){
					sprintf(input_params->basisname,"%s",line);
				}
				if(inp_indx==2){
					sscanf(line,"%lf",input_params->rxn_coord);
				}
				if(inp_indx==3){
					sprintf(input_params->method,"%s",line);
				}
				if(inp_indx==4){
					sscanf(line,"%d",input_params->roots);
				}
				inp_indx++;	
			}
		}
	}
	return input_params;
	
};

//initializes and reads
int initialize_constants(struct mol_constant* mol_constants){
	mol_constants->nmo = new int;
	mol_constants->alphae = new int;
	mol_constants->betae = new int;
	mol_constants->eSCF= new double;
	mol_constants->nuc_rep_energy = new double;
	char input_name[200];
	sprintf(input_name,"/home/settje/RFCI/molecule/%s/basis/%s/elec_constants/econst_%s_%s_%lf",mol_constants->molname,mol_constants->basisname,mol_constants->molname,mol_constants->basisname,*mol_constants->rxn_coord); 
	ifstream const_file;
	const_file.open(input_name);
	if(const_file.is_open()){
		if(!(const_file >> *mol_constants->nmo))           {cout << "Line does not exist!\n"; const_file.close();return 1;}
		if(!(const_file >> *mol_constants->alphae))        {cout << "Line does not exist!\n"; const_file.close();return 1;}
		if(!(const_file >> *mol_constants->betae))         {cout << "Line does not exist!\n"; const_file.close();return 1;}
		if(!(const_file >> *mol_constants->rxn_coord))     {cout << "Line does not exist!\n"; const_file.close();return 1;}
		if(!(const_file >> *mol_constants->nuc_rep_energy)){cout << "Line does not exist!\n"; const_file.close();return 1;}
		if(!(const_file >> *mol_constants->eSCF))          {cout << "Line does not exist!\n"; const_file.close();return 1;}
	}
	else{
		cout << "File does not exist: " << input_name << "\n";
	}

	mol_constants->astringcount = nchoosek(*mol_constants->nmo,*mol_constants->alphae);
	mol_constants->bstringcount = nchoosek(*mol_constants->nmo,*mol_constants->alphae);
	const_file.close();
	return 0;
}

//read MO integrals from file
int read_OEI(struct mol_constant * mol_constants){
	int i,j;
	mol_constants->mo_OEI = block_matrix(*mol_constants->nmo,*mol_constants->nmo);
	char input_name[200];
	sprintf(input_name,"/home/settje/RFCI/molecule/%s/basis/%s/integrals/OEI_%s_%s_%lf",mol_constants->molname,mol_constants->basisname,mol_constants->molname,mol_constants->basisname,*mol_constants->rxn_coord); 
	ifstream OEI_file;
	OEI_file.open(input_name);
	int nmo = *mol_constants->nmo;
	double buff;
	i=0; j=0;
	if(OEI_file.is_open()){
		while(OEI_file >> buff){
			mol_constants->mo_OEI[i][j] = buff;
			j++;
			if(j == nmo){
				i++;
				j = 0;
			}
			if(i == nmo){
				j=nmo;
				break;
			}
		}
	}
	else{
		cout << "File does not exist: " << input_name << "\n";
		OEI_file.close();
	}
	if(i<nmo-1 || j<nmo-1){
		cout << "Incomplete OEI file: " << input_name << "\n";
		printf("rows = %d\n cols = %d\n",i,j);
		OEI_file.close();
		return 1;
	}

	OEI_file.close();
	return 0;

}

//read MO integrals from file
int read_TEI(struct mol_constant * mol_constants){
	int i,j,k,l;
	int nmo = *mol_constants->nmo;
	mol_constants->mo_TEI = init_array(nmo*nmo*nmo*nmo);
	char input_name[200];
	sprintf(input_name,"/home/settje/RFCI/molecule/%s/basis/%s/integrals/TEI_%s_%s_%lf",mol_constants->molname,mol_constants->basisname,mol_constants->molname,mol_constants->basisname,*mol_constants->rxn_coord); 
	ifstream TEI_file;
	TEI_file.open(input_name);
	double buff;
	i=0; j=0; k=0; l=0;
	if(TEI_file.is_open()){
		while(TEI_file >> buff){
			mol_constants->mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+k*nmo+l] = buff;
			l++;
			if(l == nmo){
				k++;
				l = 0;
			}
			if(k == nmo){
				j++;
				k = 0;
			}
			if(j == nmo){
				i++;
				j = 0;
			}
			if(i == nmo){
				l=nmo;k=nmo;j=nmo;
				break;
			}
		}
	}
	else{
		cout << "File does not exist: " << input_name << "\n";
	}
	if(i<nmo-1 || j<nmo-1 || k <nmo-1 || l<nmo-1){
		cout << "Incomplete TEI file: " << input_name << "\n";
		return 1;
	}
	TEI_file.close();
	return 0;

}

//compute and store transformed OEI for use in sigma algorithm
void precompute_OEIprime(struct mol_constant * mol){
	int nmo = *mol->nmo;
	mol->mo_OEIprime = block_matrix(*mol->nmo,*mol->nmo);
	for(int k=0;k<*mol->nmo;k++){
		for(int l=0;l<*mol->nmo;l++){
			mol->mo_OEIprime[k][l]=mol->mo_OEI[k][l];
			for(int j=0;j<*mol->nmo;j++){
				mol->mo_OEIprime[k][l]-=0.5*mol->mo_TEI[k*nmo*nmo*nmo+j*nmo*nmo+j*nmo+l];

			}
		}
	}

}

void mol_const_free(struct mol_constant* mc){
	free(mc->molname);
	free(mc->basisname);
	free(mc->method);
	free(mc->alphae);
	free(mc->betae);
	free(mc->nmo);
	free(mc->roots);
	free(mc->rxn_coord);
	free(mc->eSCF);
	free(mc->nuc_rep_energy);
	free(mc->mo_TEI);
	free_block(mc->mo_OEI);
	free_block(mc->mo_OEIprime);
}

void print_integrals(FILE *output, struct mol_constant* mc){
	print_mat(mc->mo_OEI,*mc->nmo,*mc->nmo,output);
	int i,j,k,l;
	int nmo = *mc->nmo;
	for(i=0;i<nmo;i++){
		for(j=0;j<nmo;j++){
			for(k=0;k<nmo;k++){
				for(l=0;l<nmo;l++){
					fprintf(output,"%20.16lf\n",mc->mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+k*nmo+l]);
				}
			}
		}
	}
}

void print_mol_constants(struct mol_constant* mc){
	printf("molecule: %s\n",mc->molname);
	printf("basis: %s\n",mc->basisname);
	printf("nmo: %d\n",*mc->nmo);
	
}
#endif
