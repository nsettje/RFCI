/*!
** \file
** \brief This file includes the integer versions of several psi routines
** for handling arrays and matrices of doubles 
**
** David Sherrill, 1996
**
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <strings.h>

int * init_int_array(int size);
/*!
** init_int_array(): Allocates memory for one-D array of ints of dimension 
** 'size' and returns pointer to 1st element.  Zeroes all elements.
**
** Just modified the init_array() routine to do int's instead.
** This will avoid the temptation to allocate 5 integers by  
**    p = (int *) init_array(5/2), which is bad.             
**
** \param size = length of array to allocate
**
** Returns: pointer to new array
**
** C. David Sherrill
** \ingroup CIOMR
*/
int * init_int_array(int size)
{
  int *array;

  if ((array = (int *) malloc(sizeof(int)*size))==NULL) {
    fprintf(stderr,"init_array:  trouble allocating memory \n");
    fprintf(stderr,"size = %d\n",size);
    printf("UNABLE TO ALLOCATE INTEGER ARRAY. UNEDFINED BEHAVIOR ON THE WAY!\n");
    return(array);
}
  bzero(array,sizeof(int)*size);
  return(array);
};

