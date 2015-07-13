//SlaterDeterminant class taken from PSI4 for use in building the entire Hamiltonian for debugging purposes

#ifndef RFCI_slater
#define RFCI_slater
#include <stdio.h>
#include <stdlib.h>
class SlaterDeterminant {

   protected:
      unsigned nalp;
      unsigned nbet;
      unsigned char *Occs[2];

   public:
      SlaterDeterminant() { nalp=0; nbet=0; Occs[0]=NULL; Occs[1]=NULL; }
      ~SlaterDeterminant() { 
         if (Occs[0] != NULL) free(Occs[0]);
         if (Occs[1] != NULL) free(Occs[1]);
         }
      void set(unsigned int nalp, unsigned char *alpoccs, 
         unsigned int nbet, unsigned char *betoccs);
      void print(void);
      void print(FILE *outfile);
      void print_config(FILE *outfile) ;
      SlaterDeterminant& operator=(const SlaterDeterminant& s) ;
      friend int operator==(SlaterDeterminant& s1, SlaterDeterminant& s2) ;
      friend double matrix_element(SlaterDeterminant* I, SlaterDeterminant* J, double **Smo, double *pqrsINT, int nmo);
};
double matrix_element(SlaterDeterminant* I, SlaterDeterminant* J, double **mo_OEI, double *mo_TEI, int nmo);
double get_onel(int i, int j, double **mo_OEI);
double get_twoel(int i, int j, int k, int l,double *mo_TEI, int nmo);
int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J, 
   int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same, 
   int extended);
void common_orbs(int *same_alpha, int *same_beta, int cnt_alpha,
   int cnt_beta, int *common_docc, int *common_alpha_socc, 
   int *common_beta_socc, int *cnt_docc, int *cnt_alpha_socc, 
   int *cnt_beta_socc);


#endif
