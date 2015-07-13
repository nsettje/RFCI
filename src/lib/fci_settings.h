#ifndef RFCI_fci_settings
#define RFCI_fci_settings
#include <stdio.h>
#include <string.h>

struct fci_settings{
	double *FCI_CUTOFF;
	int *FCI_MAX_ITERS, *FCI_DAVIDSON_COLLAPSE;
};

struct fci_settings *read_fci_settings(){
	struct fci_settings * fci_set = (fci_settings*) malloc(sizeof(struct fci_settings));
	fci_set->FCI_CUTOFF = new double;
	fci_set->FCI_MAX_ITERS = new int;
	fci_set->FCI_DAVIDSON_COLLAPSE = new int;
	char fci_setting_file[] = "fci_settings.txt";
	FILE *settings = fopen(fci_setting_file,"r");
	if(settings == NULL){
		printf("No such file: %s\n",fci_setting_file);
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
					sscanf(line,"%lf",fci_set->FCI_CUTOFF);
				}
				if(inp_indx==1){
					sscanf(line,"%d",fci_set->FCI_MAX_ITERS);
				}
				if(inp_indx==2){
					sscanf(line,"%d",fci_set->FCI_DAVIDSON_COLLAPSE);
				}
				inp_indx++;	
			}
		}
	}

	fclose(settings);
	return fci_set;
}

void free_fci_settings(fci_settings *fs){
	free(fs->FCI_CUTOFF);
	free(fs->FCI_MAX_ITERS);
	free(fs->FCI_DAVIDSON_COLLAPSE);
}

#endif
