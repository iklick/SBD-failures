#ifndef __SBD_H__
#define __SBD_H__

/**
 * @defgroup SBD Simultaneous Block Diagonalization
 * @brief 
 * @details
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Include headers for libraries
#include "igraph.h"
#include "traces.h"


/**
 * @ingroup SBD
 * @brief Store the matrices to SBD, Ak, and the block diagonal matrices
 * @details
 */
typedef struct {
  int N; // Dimension of each matrix
  int M; // Number of matrices
  igraph_matrix_t *A; // The matrices
  igraph_matrix_t *B; // The block diagonal matrices
  double *work;
} SBD_t;

int SBD_init (SBD_t *SBD, int N, int M);
int SBD_from_graph (SBD_t *SBD, igraph_t *G, bool orbits, bool laplacian);
int SBD_free (SBD_t *SBD);

int SBD_print (SBD_t *SBD);


int SBD_lanczos (SBD_t *SBD, double epsilon, igraph_matrix_t *P);


int SBD_product (SBD_t *SBD, double *X, double *Y);
int SBD_residual (SBD_t *SBD, igraph_matrix_t *U, double *error);
int SBD_block_structure (SBD_t* SBD, igraph_vector_int_t *blocks);

/**
 * @ingroup qgraph
 * @brief The quotient graph structure stores populations, loop weights, edge weights
 * @details
 */
typedef struct {
  int nnode;    /**< Number of quotient nodes */
  int nedge;    /**< Number of quotient edges */
  int *loop;    /**< Intra-orbit degrees */
  int *pop;     /**< Orbit populations */
  int *edge;    /**< Edge array */
  int *weight;  /**< Inter-orbit degrees */
} SBD_qgraph_t;

int SBD_qgraph_init (SBD_qgraph_t *Q, int n, int m);
int SBD_qgraph_from_igraph (SBD_qgraph_t *Q, igraph_t *G);
int SBD_qgraph_free (SBD_qgraph_t *Q);

int SBD_qgraph_complete (SBD_qgraph_t *Q);

int SBD_qgraph_wire (SBD_qgraph_t *Q, igraph_t *G);

// Interface with nauty
int SBD_nauty_solve (igraph_t *G, igraph_vector_int_t *orbits);

// MBC
int SBD_belykh (igraph_t *G, igraph_vector_int_t *coloring);

// BLAS ROUTINES
// dot product
double ddot_ (int *n, double *x, int *incx, double *y, int *incy);

// vector scale
void dscal_ (int *n, double *a, double *x, int *incx);

// 2-norm
double dnrm2_ (int *n, double *x, int *incx);

// a*x + y
void daxpy_ (int *n, double *a, double *x, int *incx, double *y, int *incy);

// General matrix-vector product
void dgemv_ (char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

// Symmetric matrix-vector product
void dsymv_ (char *uplo, int *n, double *alpha, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

// matrix-matrix product update
void dsyrk_ (char *uplo, char *trans, int *n, int *m, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);

// general matrix-matrix product
void dgemm_ (char *transa, char *transb, int *n, int *m, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

// LAPACK ROUTINES
// General matrix solve
void dgesv_ (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

// General matrix eigenvalues
void dgeev_ (char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

// Symmetric matrix eigenvalues
void dsyev_ (char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

// Symmetric matrix eigenvalues with range selection
void dsyevx_ (char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *ifail, int *info);

/**
 * @ingroup arpack
 * @brief Structure to store required arrays to use Arpack
 * @details
 */
typedef struct {
  int ido;   /**< Reverse communication */
  char bmat;  /**< Whether to use generalized or normal eigenvalue problem */
  int n;   /**< Dimension of the matrix */
  char which[2];  /**< Which eigenvalues to compute */
  int nev;   /**< Number of eigenvalues to compute */
  double tol;   /**< Tolerance to accept eigenvalues */
  igraph_vector_t resid;  /**< Residual vector */
  int ncv;   /**< Dimension of the Lanczos basis */
  igraph_matrix_t v;   /**< Eigenvectors are returned here */
  igraph_vector_int_t iparam;  /**< Parameters */
  igraph_vector_int_t ipntr;   /**< Where input and ouput vectors are provided */
  igraph_vector_t workd;  /**< Workspace for vectors */
  igraph_vector_t workl;  /**< Workspace for Arpack */
  int lworkl;  /**< Length of the workspace array */
  int info;   /**< Return code */
  int rvec;  /**< whether to build the eigenvectors */
  igraph_vector_int_t select;  /**< Which eigenvectors to create */
  igraph_vector_t d;  /**< Eigenvalues are returned here */
  double sigma;  /**< Shifts (not used) */
} SBD_arpack_t;

int SBD_arpack_init (SBD_arpack_t *eig, int n, int nev);
int SBD_arpack_free (SBD_arpack_t *eig);

// ARPACK ROUTINES
void dsaupd_ (int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
void dseupd_ (int *rvec, char *howmny, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
#endif // __SBD_H__
