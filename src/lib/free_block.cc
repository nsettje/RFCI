//Routine to free a double pointer

#include <stdio.h>
void free_block(double **array){
	if(array==NULL) return;
	delete [] array[0];
	delete [] array;
}
