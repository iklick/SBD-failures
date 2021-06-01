
#include "SBD.h"

/**
 * @ingroup SBD
 * @brief Compute the matrix-vector product S*x = y without forming S
 * @param[in] SBD The SBD instance with A matrices loaded
 * @param[in] X The input vector
 * @param[out] Y The output vector
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_product (SBD_t *SBD, double *X, double *Y) {
  
  // Determine problem dimension
  int n = SBD->N;
  int N = n;
  int nn = n*n;
  double zero = 0.0;
  double pone = 1.0;
  double mone = -1.0;
  int one = 1;
  
  // Extract workspace from SBD
  double *work1 = &SBD->work[0];
  double *work2 = &SBD->work[n*n];
  
  // Initialize the Y vector to be all zeros
  memset (Y, 0, nn*sizeof(double));
  
  // Compute each matrix-matrix product in the commutators
  for (int k = 0; k < SBD->M; k++) {
    dgemm_ ("N", "N", &N, &N, &N, &pone, VECTOR(SBD->A[k].data), &N,
            X, &N, &zero, work1, &N); 
    dgemm_ ("N", "N", &N, &N, &N, &mone, X, &N, 
            VECTOR(SBD->A[k].data), &N, &pone, work1, &N);
    
    dgemm_ ("T", "N", &N, &N, &N, &pone, VECTOR(SBD->A[k].data), &N,
            work1, &N, &zero, work2, &N); 
    dgemm_ ("N", "T", &N, &N, &N, &mone, work1, &N, 
            VECTOR(SBD->A[k].data), &N, &pone, work2, &N); 
            
    daxpy_ (&nn, &pone, work2, &one, Y, &one);
  }

  return 0;
}

/**
 * @ingroup SBD
 * @brief Compute the commutator residual R <- U*A - A*U
 * @param[in] SBD The SBD instance with A matrices loaded
 * @param[in] u The matrix to check
 * @param[out] error The sum of the residual errors
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_residual (SBD_t *SBD, igraph_matrix_t *u, double *error) {
  
  // Problem dimensions
  int nn = SBD->N*SBD->N;
  
  // Constants for calling Fortran routines
  double pone = 1.0;
  double mone = -1.0;
  double zero = 0.0;
  int one = 1;
  
  // Extract some workspace
  double *work = &SBD->work[nn];
  
  // Initialize the error to zero
  *error = 0.0;
  
  // Compute the commutator
  for (int k = 0; k < SBD->M; k++) {
    dgemm_ ("N", "N", &SBD->N, &SBD->N, &SBD->N, &pone, &MATRIX(SBD->A[k],0,0), &SBD->N, &MATRIX(*u,0,0), &SBD->N, &zero, work, &SBD->N); 
    dgemm_ ("N", "N", &SBD->N, &SBD->N, &SBD->N, &mone, &MATRIX(*u,0,0), &SBD->N, &MATRIX(SBD->A[k],0,0), &SBD->N, &pone, work, &SBD->N); 
    
    // Update the error as the 2-norm of the residual matrix
    *error += dnrm2_ (&nn, work, &one);
  }

  return 0;
}

/**
 * @ingroup SBD
 * @param[in] SBD The SBD instance with block diagonal B matrices found
 * @param[out] block The block structure determined where if two elements are in the same block, they will have the same index
 * @return
 * @details
 * To determine the block structure, a graph is created which represents the set of B matrices such that two nodes, i and j, have an edge between them if at least one B matrix, say the k'th one, has a non-zero entry in B_{i,j;k}. 
 * @author Isaac Klickstein
 */
int SBD_block_structure (SBD_t *SBD, igraph_vector_int_t *block) {
  
  // Problem size
  int n = SBD->N;
  
  // Prepare the block array
  igraph_vector_int_resize (block, n);
  igraph_vector_int_fill (block, -1);
  
  // Tolerance to accept a non-zero value
  double epsilon = 1.0E-10;
  
  // The block structure is determined from a graph (not efficient)
  igraph_t graph;
  
  // Count the number of edges that will appear in the graph
  int nedge = 0;
  for (int p = 0; p < SBD->M; p++) {
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        if (fabs(MATRIX(SBD->B[p],i,j)) > epsilon) {
          nedge += 1;
        }
      }
    }
  }
  
  // Create the array of edges on the second pass
  igraph_vector_t edge;
  igraph_vector_init (&edge, 2*nedge);
  int count = 0;
  for (int p = 0; p < SBD->M; p++) {
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        if (fabs(MATRIX(SBD->B[p],i,j)) > epsilon) {
          VECTOR(edge)[2*count] = i;
          VECTOR(edge)[2*count+1] = j;
          count += 1;
        }
      }
    }
  }
  
  // Create the graph and remove multi-edges
  igraph_create (&graph, &edge, n, IGRAPH_UNDIRECTED);
  igraph_simplify (&graph, 1, 1, NULL);
  
  // Determine the connected components of the graph
  count = 0;
  for (int i = 0; i < n; i++) {
    if (VECTOR(*block)[i] == -1) {
      VECTOR(*block)[i] = count;
      igraph_subcomponent (&graph, &edge, i, IGRAPH_ALL);
      for (int j = 0; j < igraph_vector_size (&edge); j++) {
        VECTOR(*block)[(int)VECTOR(edge)[j]] = count;
      }
      count += 1;
    }
  }
  
  // Free the local memory
  igraph_destroy (&graph);
  igraph_vector_destroy (&edge);

  return 0;
}
