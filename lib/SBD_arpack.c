
#include "SBD.h"

/**
 * @ingroup arpack
 * @brief Initialize an arpack interface
 * @param[out] eig The arpack interface to initialize
 * @param[in] n The dimension of the matrix 
 * @param[in] nev The number of eigenvalues to compute (must be less than n)
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_arpack_init (SBD_arpack_t *eig, int n, int nev) {

  eig->ido = 0;
  eig->bmat = 'I';
  eig->n = n*n;
  
  eig->which[0] = 'S';
  eig->which[1] = 'M';
  
  eig->nev = nev;
  eig->tol = 1.0E-14;
  
  igraph_vector_init (&eig->resid, n*n);
  eig->ncv = nev+10;
  
  igraph_matrix_init (&eig->v, n*n, eig->ncv);
  igraph_vector_int_init (&eig->iparam, 11);
  igraph_vector_int_init (&eig->ipntr, 14);
  
  VECTOR(eig->iparam)[0] = 1;
  VECTOR(eig->iparam)[1] = 0;
  VECTOR(eig->iparam)[2] = n*n*n;
  VECTOR(eig->iparam)[3] = 1;
  VECTOR(eig->iparam)[4] = 0;
  VECTOR(eig->iparam)[5] = 0;
  VECTOR(eig->iparam)[6] = 1;
  
  igraph_vector_init (&eig->workd, 3*n*n);
  igraph_vector_init (&eig->workl, eig->ncv * eig->ncv + 8*eig->ncv + 10);
  eig->lworkl = eig->ncv * eig->ncv + 8 * eig->ncv + 10;
  
  igraph_vector_int_init (&eig->select, eig->ncv); 
  igraph_vector_init (&eig->d, nev);
  
  eig->sigma = 0;
  eig->rvec = 1;
  
  eig->info = 0;
  
  return 0;
}

/**
 * @ingroup arpack
 * @brief Free memory allocated for an arpack interface
 * @param[in] eig The arpack interface to free
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_arpack_free (SBD_arpack_t *eig) {
  igraph_vector_destroy (&eig->resid);
  igraph_matrix_destroy (&eig->v);
  igraph_vector_int_destroy (&eig->iparam);
  igraph_vector_int_destroy (&eig->ipntr);
  igraph_vector_destroy (&eig->workd);
  igraph_vector_destroy (&eig->workl);
  igraph_vector_int_destroy (&eig->select);
  igraph_vector_destroy (&eig->d);
  return 0;
}
