/* This file defines a structure that carries relevant settings for the RFCI routine as read from the file rfci_settings.txt in ~/RFCI.
 *
 * output:
 * struct of RFCI settings to be read by other routines
 *
 * input:
 * rfci_settings.txt file
 */

#ifndef RRFCI_rrfci_settings
#define RRFCI_rrfci_settings
#include <string.h>
struct rfci_settings{
	double *FCI_CUTOFF,*RFCI_CUTOFF;
	int *RFCI_MAX_WFN_ITERS, *RFCI_MAX_DAV_ITERS, *RFCI_DAVIDSON_COLLAPSE, *FCI_DAVIDSON_COLLAPSE;
};
struct rfci_settings *read_rfci_settings(){
	struct rfci_settings * rfci_set = (rfci_settings*) malloc(sizeof(struct rfci_settings));
	rfci_set->RFCI_CUTOFF = new double; //Davidson convergence criterion
	rfci_set->RFCI_MAX_WFN_ITERS = new int; //Maximum number of terms to include in wfn
	rfci_set->RFCI_MAX_DAV_ITERS = new int; //Maximum # of Davidson iterations allowed per diagonalization
	rfci_set->RFCI_DAVIDSON_COLLAPSE = new int; //Whether or not to allow subspace collapse in Davidson
	char rfci_setting_file[] = "rfci_settings.txt";
	FILE *settings = fopen(rfci_setting_file,"r");
	if(settings == NULL){
		printf("No such file: %s\n",rfci_setting_file);
	}
	if(settings!= NULL){
		char line[128];
		int inp_indx = 0;
		while(fgets(line,sizeof line,settings)){
			if(line[0] != '#'){
				char *pos;
				if((pos=strchr(line,'\n')) != NULL){ *pos = '\0';}
				sprintf(line,"%s",line); //strips the '\n' off the end of the line
				if(inp_indx==0){
					sscanf(line,"%lf",rfci_set->RFCI_CUTOFF);
				}
				if(inp_indx==1){
					sscanf(line,"%d",rfci_set->RFCI_MAX_WFN_ITERS);
				}
				if(inp_indx==2){
					sscanf(line,"%d",rfci_set->RFCI_MAX_DAV_ITERS);
				}
				if(inp_indx==3){
					sscanf(line,"%d",rfci_set->RFCI_DAVIDSON_COLLAPSE);
				}
				inp_indx++;	
			}
		}
	}

	fclose(settings);
	return rfci_set;
}

void free_rfci_settings(rfci_settings *fs){
	free(fs->RFCI_CUTOFF);
	free(fs->RFCI_MAX_WFN_ITERS);
	free(fs->RFCI_MAX_DAV_ITERS);
	free(fs->RFCI_DAVIDSON_COLLAPSE);
}
#endif
