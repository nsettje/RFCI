/* This is a function that allocates a 2D array.
 *
 * output: 
 * double **matrix of size n-by-m
 *
 * input: 
 * n (rows)
 * m (cols)
 * */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>
#include <unistd.h>
//#include "memory.h"
#ifdef _POSIX_MEMLOCK
#include <sys/mman.h>
#endif

#ifdef HAVE_MKL
#ifdef HAVE_MKL_MALLOC

extern "C" {
void* MKL_malloc(size_t size, int align);
void MKL_free(void *ptr);
}

#endif
#endif


/*! Taken from PSI4
** block_matrix(): Allocate a 2D array of doubles using contiguous memory
**
** Allocates a contiguous block of memory for an array of
** doubles, allocates an array of pointers to the beginning of each row and
** returns the pointer to the first row pointer.  This allows transparent
** 2d-array style access, but keeps memory together such that the matrix
** could be used in conjunction with FORTRAN matrix routines.
**
** Allocates memory for an n x m matrix and returns a pointer to the
** first row.
**
** \param n = number of rows (unsigned long to allow large matrices)
** \param m = number of columns (unsigned long to allow large matrices)
** \param memlock = optional bool indicating whether to lock memory
**   into physical RAM or not, and available only where _POSIX_MEMLOCK
**   is defined. Defaults to false if not specified.
**
** Returns: double star pointer to newly allocated matrix
**
** T. Daniel Crawford
** Sometime in 1994
**
** Based on init_matrix() from libciomr
** \ingroup CIOMR
*/

double ** block_matrix(int n, int m)
{
	bool memlock = true;
    double **A=NULL;
    double *B=NULL;
    unsigned long int i;

    if(!m || !n) return(static_cast<double **>(0));

    A = new double*[n];
    if (A==NULL) {
        fprintf(stderr,"block_matrix: trouble allocating memory \n");
        fprintf(stderr,"n = %d\n",n);
    }

    B = new double[n*m];
    if (B == NULL) {
        fprintf(stderr,"block_matrix: trouble allocating memory \n");
        fprintf(stderr,"m = %d\n",m);
    }
    memset(static_cast<void*>(B), 0, m*n*sizeof(double));

    for (i = 0; i < n; i++) {
        A[i] = &(B[i*m]);
    }

#ifdef _POSIX_MEMLOCK
    if (memlock) {

        char* addr = (char*) B;
        unsigned long size = m*n*(unsigned long)sizeof(double);
        unsigned long page_offset, page_size;

        page_size = sysconf(_SC_PAGESIZE);
        page_offset = (unsigned long) addr % page_size;

        addr -= page_offset;  /* Adjust addr to page boundary */
        size += page_offset;  /* Adjust size with page_offset */

        if ( mlock(addr, size) ) {  /* Lock the memory */
            //fprintf(stderr,"block_matrix: trouble locking memory \n");
            fflush(stderr);
        }

        addr = (char*) A;
        size = n*(unsigned long)sizeof(double*);

        page_offset = (unsigned long) addr % page_size;

        addr -= page_offset;  /* Adjust addr to page boundary */
        size += page_offset;  /* Adjust size with page_offset */

        if ( mlock(addr, size) ) {  /* Lock the memory */
            //fprintf(stderr,"block_matrix: trouble locking memory \n");
            fflush(stderr);
        }
    }
#endif

    return(A);
}


/*!
** free_block(): Free a block matrix
**
** \param array = pointer to matrix to be freed
**
** Returns: none
**
** \ingroup CIOMR
*/

