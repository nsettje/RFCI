//

/*
** \file
** \brief Initialize an array of doubles
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <strings.h>


/*!
** init_array(): This function initializes an array of doubles of
** length 'size' and returns a pointer to the first element
**
** \param size = length of array (unsigned long to allow large arrays)
**
** Returns: pointer to new array
**
** \ingroup CIOMR
*/
double * init_array(int size)
{
  double *array;

  if ((array = (double *) malloc(size*(unsigned long int)sizeof(double)))
    == NULL) {
    fprintf(stderr,"init_array: trouble allocating memory \n");
    //fprintf(stderr,"size = %ld\n",size);
  }
  bzero(array,size*(unsigned long int)sizeof(double));
  return(array);
}


