
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"
#include "mmio.c"
#include <cilk/cilk.h>
struct timeval start, end;
double ex_time;
//function coo2csc converts COO format to CSC

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}



int main(int argc, char *argv[]) {

    /***********************************************/   
    /*Start reading Matrix Market files(COO format)*/
    /***********************************************/

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t m, n, nz;
    uint32_t i,j,k,l, *coo_row, *coo_col;
    double *val;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */
    printf("%d \n",n);
    coo_row = (uint32_t *)malloc(nz * sizeof(uint32_t));
    coo_col = (uint32_t *)malloc(nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode))
    {
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &coo_row[i], &coo_col[i], &val[i]);
        coo_row[i]--;  /* adjust from 1-based to 0-based */
        coo_col[i]--;
    }
    }
    else
    {
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d\n", &coo_row[i], &coo_col[i]);
        val[i]=1;
        coo_row[i]--;  /* adjust from 1-based to 0-based */
        coo_col[i]--;
    }
    }



    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, m, n, nz);
    
    /************************************************/   
    /*Finish reading Matrix Market files(COO format)*/
    /************************************************/
 



    uint32_t * csc_row = (uint32_t *)malloc(nz     * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc((n + 1) * sizeof(uint32_t));

    uint32_t isOneBased = 0;

    // Call coo2csc for isOneBase false
    coo2csc(csc_row, csc_col,
          coo_row, coo_col,
          nz, n,
          isOneBased);

    // Verify output

    /*for (i = 0; i < n + 1 ; i++) {
    printf("%d ", csc_col[i]);
    }
    printf("\n");
    for (uint32_t i = 0; i < nz; i++) {
    	printf("%d ", csc_row[i]);
    }*/
    printf("\n");
    /*Initializing of c3 vector(number of triangles that a node belongs to*/
    uint32_t * c3 = (uint32_t *)malloc(n * sizeof(uint32_t));
    for (uint32_t i=0;i<n;i++){
		c3[i]=0;
    }
    uint32_t fir,mid;
    gettimeofday(&start,NULL); //Start timing the computation
    /********************************************/   
    /*Now parallel with CILK V3 algorithm starts*/
    /********************************************/


    cilk_for (i=0;i<n+1;i++){
    	 for(j=csc_col[i];j<csc_col[i+1];j++){
		for(k=csc_col[csc_row[j]];k<csc_col[csc_row[j]+1];k++){
			/*Binary search to check if k the neighbour of j is also a neighbour for i*/	
                			fir = j+1;
                			l = csc_col[i+1]-1;
                			mid = (fir+l)/2;
                			while (fir <= l) {
                  			if (csc_row[mid] == csc_row[k]) {
                      				#pragma omp critical
		       			c3[i]++;
		       			c3[csc_row[mid]]++;
                       			c3[csc_row[j]]++;
                      				break;
					}
                  			else if (csc_row[mid] < csc_row[k]){
                    				fir = mid + 1;
                  			}        
                  			else{
                   				l = mid - 1;
                  			    }		
					mid = (fir+l)/2;
                			}	
		}
	}
    }

    /**********************************************/   
    /*Now parallel with CILK V3 algorithm finishes*/
    /**********************************************/

    gettimeofday(&end,NULL); //Stop timing the computation    
    ex_time=(end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);
    printf( "V3 took %lf seconds to execute.\n", ex_time); 
    /*Print the number of triangles for each node*/
    for(i=0;i<n;i++){
	printf("%d %d \n",i,c3[i]);
    }

 




    /* cleanup variables */
    free( csc_row );
    free( csc_col );
    free( coo_col );
    free( coo_row );
    free( c3 );

    return 0;
}
