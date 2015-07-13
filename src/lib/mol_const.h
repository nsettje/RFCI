//Actual definition of mol_constant struct
//In a separate header file to avoid clashing symbols when including this in many files

#ifndef RFCI_mol_const
#define RFCI_mol_const
struct mol_constant {
	char * molname;
	char *basisname;
	char *method; 
	double *rxn_coord, *nuc_rep_energy, *eSCF;	
	int *nmo, *alphae, *betae, *roots;
	int astringcount, bstringcount;
	double **mo_OEI;
	double **mo_OEIprime;
	double *mo_TEI;

}; 
#endif
