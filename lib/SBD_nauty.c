
#include "SBD.h"

/**
 * @ingroup SBD
 * @brief Determine the orbits of the automorphism group using nauty
 * @param[in] G The graph to analyze
 * @param[out] orbits The orbit labels of the nodes in G
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_nauty_solve (igraph_t *G, igraph_vector_int_t *orbits) {
  
  // Problem dimensions 
  int n = igraph_vcount (G);
  int m = igraph_ecount (G);
  
  // Local storage to construct the form of the graph used by nauty
  igraph_vector_t result;
  igraph_vector_init (&result, n);
  
  // Options to pass to Traces
  static DEFAULTOPTIONS_TRACES (topt);
  topt.verbosity = 0;
  
  // Other local memory used for Traces
  TracesStats tstat;
  sparsegraph g;
  
  // Set the graph size for the nauty version of the graph
  g.nv = n;
  g.nde = 2*m;
  
  // Allocate memory for the nauty graph
  g.d = (int*) malloc (g.nv*sizeof(int));
  g.v = (size_t*) malloc ((g.nv+1)*sizeof(size_t));
  g.e = (int*) malloc (g.nde*sizeof(int));
  g.w = NULL;
  g.vlen = g.nv;
  g.elen = g.nde;
  g.dlen = g.nv;
  g.wlen = 0;
  
  // Loop over nodes
  //   The graph is stored in compressed sparse format
  g.v[0] = 0;
  for (int i = 0; i < n; i++) {
    igraph_neighbors (G, &result, i, IGRAPH_ALL);
    g.d[i] = igraph_vector_size (&result);
    g.v[i+1] = g.v[i] + g.d[i];
    for (int j = 0; j < g.d[i]; j++) {
      g.e[g.v[i]+j] = (int)VECTOR(result)[j];
    }
  }
  
  // More local memory for Traces
  int *lab = (int*) calloc (g.nv, sizeof(int));
  int *ptn = (int*) calloc (g.nv, sizeof(int));
  int *orb = (int*) calloc (g.nv, sizeof(int));
  
  // Call Traces
  Traces (&g, lab, ptn, orb, &topt, &tstat, NULL);
  
  // After returning from Traces, the orbits are labeled according to the lowest indexed member node (so orbit labels are not contiguous)
  // The orbit indices are shifted to be contiguous (from 0 to p-1)
  for (int i = 0; i < n; i++) {
    lab[i] = -1;
  }
  
  int count = 0;
  for (int i = 0; i < n; i++) {
    if (lab[orb[i]] == -1) {
      lab[orb[i]] = count;
      count += 1;
    }
  }
  
  igraph_vector_int_resize (orbits, n);
  for (int i = 0; i < n; i++) {
    VECTOR(*orbits)[i] = lab[orb[i]];
  }
  
  
  // Free the local memory
  free (g.d);
  free (g.v);
  free (g.e);
  igraph_vector_destroy (&result);
  free (lab);
  free (ptn);
  free (orb);
  return 0;
}

/**
 * @ingroup SBD
 * @brief Compute the minimum balanced coloring of a graph
 * @param[in] G The graph 
 * @param[out] coloring The coloring computed
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_belykh (igraph_t *G, igraph_vector_int_t *coloring) {
  int N = igraph_vcount (G);
  int M = igraph_ecount (G);
  igraph_vector_int_resize (coloring, N);
  igraph_vector_int_fill (coloring, -1);
  
  int v1, v2;
  
  // LOCAL MEMORY
  igraph_matrix_int_t D;
  igraph_vector_int_t c;
  
  igraph_vector_int_init (&c, N);
  igraph_vector_int_fill (&c, 0);
  igraph_matrix_int_init (&D, N, N);
  
  int ncolor = 1;
  
  int flag = 1;
  while (flag) {
    
    // Initialize the coloring to create array
    igraph_vector_int_fill (coloring, -1);
  
    // Create the cluster degree matrix
    igraph_matrix_int_fill (&D, 0);
    for (int k = 0; k < M; k++) {
      igraph_edge (G, k, &v1, &v2);
      MATRIX(D,v1,VECTOR(c)[v2]) += 1;
      MATRIX(D,v2,VECTOR(c)[v1]) += 1;
    }
    
    // Perform the refinement
    int color = 0;
    for (int i = 0; i < N; i++) {
      if (VECTOR(*coloring)[i] == -1) {
        VECTOR(*coloring)[i] = color;
        for (int j = i+1; j < N; j++) {
          if (VECTOR(c)[j] == VECTOR(c)[i]) {
            int eq = 1;
            for (int k = 0; k < ncolor; k++) {
              if (MATRIX(D,i,k) != MATRIX(D,j,k)) {
                eq = 0;
                break;
              }
            }
            if (eq) {
              VECTOR(*coloring)[j] = color;
            }
          }
        }
        color += 1;
      }
    }
    
    ncolor = color;
    
    // Check if coloring and c are equal, i.e., no changes
    if (igraph_vector_int_all_e (coloring, &c)) flag = 0;
    else igraph_vector_int_update (&c, coloring);
  }
  
  return 0;
}
