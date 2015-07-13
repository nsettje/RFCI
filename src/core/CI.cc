#include <stdio.h>
#include "lib/IO.h"
#include "lib/FCI.h"
#include "lib/RFCI.h"
#include "lib/mkl_wrapper.h"
#include "lib/mkl_test.h"

int main(int argc, char* argv[]){
//BEGIN INITIALIZATIONS
	//parse inputs: CI [input_file.txt] [output_file.txt] with defaults "input.txt" and "output.txt"
	char *input;
	char *output;
	if(argc == 1){
		input = "input.txt";
		output = "output.txt";
	}
	else if(argc == 2){
		input = argv[1];
		output = "output.txt";
	}
	else if(argc==3){
		input = argv[1];
		output = argv[2];
	}
	//initialize constants and integrals in struct
	struct mol_constant *mol = mem_input(input); //read input
	initialize_constants(mol); //read econst file
	read_OEI(mol); //read OEI file
	read_TEI(mol);//read TEI file
	//pre-compute intermediate quantities to use in the sigma builds
	precompute_OEIprime(mol);
	print_mol_constants(mol);
//END INITIALIZATIONS

	FILE * outfile = fopen(output,"w");

//BEGIN METHOD
	//FCI method
	if(strcmp(mol->method,"FCI")==0){
		fci_davidson(mol,output,1);
	}

	if(strcmp(mol->method,"RFCI")==0){
		rfci_slow_method(mol,output,1);
	}
	

	mol_const_free(mol);
//END METHOD
	//write outpu
	fclose(outfile);
	//clean up memory



	return 0;
}
