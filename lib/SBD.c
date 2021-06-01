
#include "SBD.h"

/**
 * @ingroup SBD
 * @brief Initialize an SBD instance
 * @param[out] SBD The SBD instance to initialize
 * @param[in] N The dimension of the matrices to SBD
 * @param[in] M The number of matrices to SBD
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_init (SBD_t *SBD, int N, int M) {
  
  // Check to see if N and M are positive
  if (N <= 0 || M <= 0) {
    printf("SBD requires positive dimensions\n");
    return 1;
  }
  
  // Set the dimensions
  SBD->N = N;
  SBD->M = M;
  
  // Allocate memory for the matrix arrays
  SBD->A = (igraph_matrix_t*) calloc (SBD->M, sizeof(igraph_matrix_t));
  SBD->B = (igraph_matrix_t*) calloc (SBD->M, sizeof(igraph_matrix_t));
  for (int k = 0; k < M; k++) {
    igraph_matrix_init (&SBD->A[k], N, N);
    igraph_matrix_init (&SBD->B[k], N, N);
  }
  
  // Create some workspace for other routines
  SBD->work = (double*) calloc (4*SBD->N*SBD->N + 2*SBD->M*SBD->N, sizeof(double));
  
  return 0;
}

/**
 * @ingroup SBD
 * @brief Initialize an SBD instance as a graph adjacency matrix and its orbit indicator matrices.
 * @param[out] SBD The SBD instance to initialize
 * @param[in] G The graph used to create the SBD instance
 * @param[in] orb Whether to use the orbits of the automorphism graph or use the minimal balanced coloring
 * @param[in] laplacian If true use the Laplacian, if false use the adjacency matrix
 * @return 
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_from_graph (SBD_t *SBD, igraph_t *G, bool orb, bool laplacian) {
  
  // Determine the size of the graph
  int n = igraph_vcount (G);
  
  // Initialize storage for the orbits
  igraph_vector_int_t orbits;
  igraph_vector_int_init (&orbits, n);
  
  // Use nauty to find the orbits of the automorphism group
  if (orb) {
    SBD_nauty_solve (G, &orbits);
  }
  else {
    SBD_belykh (G, &orbits);
  }
  // Determine the number of orbits
  int M = igraph_vector_int_max (&orbits) + 1;
  
  // Initialize the SBD with known size
  SBD_init (SBD, n, M+1);
  
  // Set the first matrix in the SBD to be the adjacency matrix of the graph
  if (laplacian) {
    igraph_laplacian (G, &SBD->A[0], NULL, 0, NULL);
  }
  else {
    igraph_get_adjacency (G, &SBD->A[0], IGRAPH_GET_ADJACENCY_BOTH, 0);
  }
  // Make sure all the other matrices are initially all zero
  for (int k = 1; k < SBD->M; k++) {
    igraph_matrix_null (&SBD->A[k]);
  }
  
  // Build the orbit indicator matrices
  for (int k = 0; k < SBD->N; k++) {
    MATRIX(SBD->A[VECTOR(orbits)[k]+1], k, k) = 1.0;
  }
  
  // Free the locally allocated memory
  igraph_vector_int_destroy (&orbits);
  return 0;
}

/** 
 * @ingroup SBD
 * @brief Free the memory allocated for an SBD instance
 * @param[in] SBD The SBD instance to free
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_free (SBD_t *SBD) {
  for (int k = 0; k < SBD->M; k++) {
    igraph_matrix_destroy (&SBD->A[k]);
    igraph_matrix_destroy (&SBD->B[k]);
  }
  free (SBD->A);
  free (SBD->B);
  free (SBD->work);
  return 0;

}

/**
 * @ingroup SBD
 * @brief Print an SBD instance
 * @param[in] SBD The SBD instance to print to stdout
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_print (SBD_t *SBD) {
  printf("SBD Matrix:\n");
  for (int k = 0; k < SBD->M; k++) {
    printf("  A[%3i] = [\n", k);
    igraph_matrix_print (&SBD->A[k]);
    printf("]\n\n");
  }
  return 0;
}
