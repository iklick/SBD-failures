
#include "SBD.h"

/**
 * @ingroup SBD
 * @brief Solve the SBD problem using the Lanczos method to determine the null space of the S matrix without explicitly constructing it.
 * @param[inout] SBD The SBD instant to solve with A matrices loaded.
 * @param[in] epsilon The threshold to accept an eigenvalue as being zero
 * @param[out] P The matrix which simultaneously block diagonalizes each of the A matrices
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_lanczos (SBD_t *SBD, double epsilon, igraph_matrix_t *P) {
  
  // Get the problem size
  int n = SBD->N;
  int nn = SBD->N*SBD->N;
  
  // Store some constants for calling Fortran routines
  int one = 1;
  double zero = 0.0;
  double pone = 1.0;
  
  // Initialize the arpack interface
  SBD_arpack_t eig;
  SBD_arpack_init (&eig, SBD->N, 10);
  eig.ido = 0;
  
  // Storage for the matrix-vector product input and output
  double *x, *y;
  
  // Begin the reverse communication loop
  int idx = 0;
  while (1) {
  
    // DSAUPD performs each Lanczos step, returning when a matrix-vector product is required
    dsaupd_ (&eig.ido, &eig.bmat, &eig.n, eig.which, &eig.nev, &eig.tol,
             VECTOR(eig.resid), &eig.ncv, VECTOR(eig.v.data), &eig.n,
             VECTOR(eig.iparam), VECTOR(eig.ipntr), VECTOR(eig.workd),
             VECTOR(eig.workl), &eig.lworkl, &eig.info);
    
    // If the reverse communication requestions a matrix-vector product, provide it         
    if (eig.ido == 1 || eig.ido == -1) {
      x = &VECTOR(eig.workd)[VECTOR(eig.ipntr)[0]-1];
      y = &VECTOR(eig.workd)[VECTOR(eig.ipntr)[1]-1];
      SBD_product (SBD, x, y);
    }
    else {
      break;
    }     
    
    //printf("  Step: %i\n", idx);
    idx += 1;
    
  }
  
  
  // Post-processing step
  dseupd_ (&eig.rvec, "A", VECTOR(eig.select), VECTOR(eig.d), 
           VECTOR(eig.v.data), &eig.n, &eig.sigma, &eig.bmat, &eig.n,
           eig.which, &eig.nev, &eig.tol, VECTOR(eig.resid), &eig.ncv,
           VECTOR(eig.v.data), &eig.n, VECTOR(eig.iparam), VECTOR(eig.ipntr),
           VECTOR(eig.workd), VECTOR(eig.workl), &eig.lworkl, &eig.info);
  
  // Create some matrices from the SBD instance workspace
  igraph_matrix_t U, A;
  igraph_matrix_view (&U, SBD->work, SBD->N, SBD->N);
  igraph_matrix_view (&A, &SBD->work[nn], SBD->N, SBD->N);
  igraph_matrix_null (&A);
  
  int m = 0;
  for (int k = 0; k < eig.nev; k++) {
    //printf("Eig[%3i] = %E\n", k, VECTOR(eig.d)[k]);
    if (VECTOR(eig.d)[k] < epsilon) {
      m += 1;
    }
  }
  
  // Create a random U matrix from the nullspace basis
  for (int k = 0; k < m; k++) {
    double c = 1.0 / (double)m;
    daxpy_ (&nn, &c, &MATRIX(eig.v,0,k), &one, &VECTOR(A.data)[0], &one);
  }
  
  
  // If U is in the nullspace of A, then so U^T + U which is symmetric and 'nicer'
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      MATRIX(U,i,j) = ( MATRIX(A,i,j) + MATRIX(A,j,i) )/2.0;
    }
  }
  
  // Compute the commutating residual error to ensure U is a good choice
  double error;
  SBD_residual (SBD, &U, &error);
  printf("Residual: %E\n", error);
  
  // Create some workspace for the dense eigenvalue solver coming up
  int lwork = 10*(n+5);
  double *work = (double*) calloc (lwork, sizeof(double));
  double *w = (double*) calloc (n, sizeof(double));
  
  int info = 0;

  // The SBD matrix P is the matrix of eigenvectors of U  
  igraph_matrix_update (P, &U);
  dsyev_ ("V", "U", &n, VECTOR(P->data), &n, w, work, &lwork, &info);

  // Check return code from DSYEV    
  if (info != 0) {
    printf("DSYEV failed with code %i\n", info);
  }
  
  // Perform the SBD and store them in the array of B matrices
  for (int k = 0; k < SBD->M; k++) {
    dgemm_ ("N", "N", &n, &n, &n, &pone, VECTOR(SBD->A[k].data), &n, VECTOR(P->data), &n, &zero, VECTOR(U.data), &n);
    dgemm_ ("T", "N", &n, &n, &n, &pone, VECTOR(P->data), &n,VECTOR(U.data) , &n, &zero, VECTOR(SBD->B[k].data), &n);
  }
  
  // Free the local memory
  SBD_arpack_free (&eig);
  free (work);
  free (w);

  return 0;
}

