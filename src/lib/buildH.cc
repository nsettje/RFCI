/*This is a routine that builds the entire full CI Hamiltonian (only for debugging purpose)
 *
 * output:
 * double **H (Hamiltonian, astringcount by bstringcount)
 *
 * input:
 * mo_OEI (MO one-electron integrals)
 * mo_TEI (MO two-electron integrals)
 * alphae (# alpha electrons)
 * betae  (#  beta electrons)
 * nmo 	  (# MOs)
 *
 * */

#include "lib/hamiltonian.h"
//Builds full-size Hamiltonian explicitly according to Slater's rules
//#include "lib/hamiltonian.h"
void build_full_Hamiltonian(double **H, double **mo_OEI, double *mo_TEI, int alphae, int betae, int nmo){
	int astringcount=nchoosek(nmo,alphae);
	int bstringcount=nchoosek(nmo,betae);
	int N = astringcount*bstringcount;
	//strings of occupied determinants for a, b electrons
	unsigned char astringR[alphae], astringC[alphae], bstringR[betae], bstringC[betae];
	SlaterDeterminant Row;
	SlaterDeterminant Col;
	for(int row=0;row<N;row++){ //loop over H rows
		for(int i=0;i<alphae;i++){
			astringR[i]=i;
		}
		for(int i=0;i<betae;i++){
			bstringR[i]=i;
		}
		int * indexR = init_int_array(2);
		//increment index up to appropriate level of excitation for a or b in row determinant
		for(int r_index=0;r_index<row;r_index++){ //indexes string row
			indexR[1]++;
			if(indexR[1]==astringcount){
				indexR[1]=0;
				indexR[0]++;
			}
		}
		//compute next combination for a and b
		for(int comb=0;comb<indexR[0];comb++){
			next_combination_char(astringR,nmo,alphae);
		}
		for(int comb=0;comb<indexR[1];comb++){
			next_combination_char(bstringR,nmo,alphae);
		}
		free(indexR);
		Row.set(alphae,astringR,betae,bstringR);
		for(int col=0;col<N;col++){ //loop over H cols
			for(int i=0;i<alphae;i++){
				astringC[i]=i;
			}
			for(int i=0;i<betae;i++){
				bstringC[i]=i;
			}
			int * indexC = init_int_array(2);
			//increment index up to appropriate level of excitation for a or b in col determinant
			for(int c_index=0;c_index<col;c_index++){ //indexes col row
				indexC[1]++;
				if(indexC[1]==astringcount){
					indexC[1]=0;
					indexC[0]++;
				}

			}
			for(int comb=0;comb<indexC[0];comb++){
				next_combination_char(astringC,nmo,alphae);
			}
			for(int comb=0;comb<indexC[1];comb++){
				next_combination_char(bstringC,nmo,betae);
			}
					free(indexC);
			Col.set(alphae,astringC,betae,bstringC);
			//compute matrix element
			//printf("%d %d\n",row,col);
			H[row][col]=matrix_element(&Row,&Col,mo_OEI,mo_TEI,nmo);
		}
	}

}

//computes a single element of H, specifically H(row,col)
double build_single_Hamiltonian_element(int row, int col, double **mo_OEI, double *mo_TEI, int alphae, int betae, int nmo){
	int astringcount=nchoosek(nmo,alphae);
	int bstringcount=nchoosek(nmo,betae);
	int N = astringcount*bstringcount;
	//strings of occupied determinants for a, b electrons
	unsigned char astringR[alphae], astringC[alphae], bstringR[betae], bstringC[betae];
	SlaterDeterminant Row;
	SlaterDeterminant Col;
	for(int i=0;i<alphae;i++){
		astringR[i]=i;
	}
	for(int i=0;i<betae;i++){
		bstringR[i]=i;
	}
	int * indexR = init_int_array(2);
	//increment index up to appropriate level of excitation for a or b in row determinant
	for(int r_index=0;r_index<row;r_index++){ //indexes string row
		indexR[1]++;
		if(indexR[1]==astringcount){
			indexR[1]=0;
			indexR[0]++;
		}
	}
	//compute next combination for a and b
	for(int comb=0;comb<indexR[0];comb++){
		next_combination_char(astringR,nmo,alphae);
	}
	for(int comb=0;comb<indexR[1];comb++){
		next_combination_char(bstringR,nmo,alphae);
	}
	free(indexR);
	Row.set(alphae,astringR,betae,bstringR);
	for(int i=0;i<alphae;i++){
		astringC[i]=i;
	}
	for(int i=0;i<betae;i++){
		bstringC[i]=i;
	}
	int * indexC = init_int_array(2);
	//increment index up to appropriate level of excitation for a or b in col determinant
	for(int c_index=0;c_index<col;c_index++){ //indexes col row
		indexC[1]++;
		if(indexC[1]==astringcount){
			indexC[1]=0;
			indexC[0]++;
		}

	}
	for(int comb=0;comb<indexC[0];comb++){
		next_combination_char(astringC,nmo,alphae);
	}
	for(int comb=0;comb<indexC[1];comb++){
		next_combination_char(bstringC,nmo,betae);
	}
			free(indexC);
	Col.set(alphae,astringC,betae,bstringC);
	//compute matrix element
	double val=matrix_element(&Row,&Col,mo_OEI,mo_TEI,nmo);
	return(val);

}

