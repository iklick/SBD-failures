
#include "SBD.h"


static int gcd (int a, int b) {
  if (a == 0) {
    return b;
  }
  return gcd(b%a,a);
}

static int random_factor (int c) {

  int limit = (int)floor(sqrt((double)c));
  
  int *factors = (int*) calloc (2*limit, sizeof(int));
  
  int count = 0;
  for (int i = 1; i <= limit; i++) {
    if (c%i == 0) {
      if (c/i == i) {
        factors[count] = i;
        count += 1;
      }
      else {
        factors[count] = i;
        factors[count+1] = c/i;
        count += 2;
      }
    }
  }
  
  int choice = rand() % count;
  
  int d = factors[choice];
  
  free (factors);
  return d;
}

/**
 * @ingroup qgraph
 * @brief Initialize a quotient graph
 * @param[out] Q The quotient graph to create
 * @param[in] n The number of quotient nodes
 * @param[in] m The number of quotient edges
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_qgraph_init (SBD_qgraph_t *Q, int n, int m) {
  
  Q->nnode = n;
  Q->nedge = m;
  Q->loop = (int*) calloc (Q->nnode, sizeof(int));
  Q->pop  = (int*) calloc (Q->nnode, sizeof(int));
  Q->edge = (int*) calloc (2*m, sizeof(int));
  Q->weight = (int*) calloc (2*m, sizeof(int));
  
  
  return 0;
}

/**
 * @ingroup qgraph
 * @brief Create the structure of a quotient graph from an igraph
 * @param[out] Q The quotient graph to initialize
 * @param[in] G The graph used for the structure of Q
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_qgraph_from_igraph (SBD_qgraph_t *Q, igraph_t *G) {
  int n = igraph_vcount (G);
  int m = igraph_ecount (G);
  
  SBD_qgraph_init (Q, n, m);
  int v, u;
  for (int i = 0; i < m; i++) {
    igraph_edge (G, i, &v, &u);
    Q->edge[2*i]   = v;
    Q->edge[2*i+1] = u;
  }
  
  return 0;
}

/**
 * @ingroup qgraph
 * @brief Free the memory allocated by a quotient graph
 * @param[in] Q The quotient graph to free
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_qgraph_free (SBD_qgraph_t *Q) {
  
  free (Q->loop);
  free (Q->pop);
  free (Q->edge);
  free (Q->weight);
  
  return 0;
}

/**
 * @ingroup qgraph
 * @brief Given a quotient graph with edges, pops, and loops set, find a set of feasible quotient edge weights
 * @param[inout] Q The quotient graph for which edge weights are determined
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_qgraph_complete (SBD_qgraph_t *Q) {

  // Variables for orbit index, population
  int o1, o2, n1, n2;
  int c, d;
  
  for (int i = 0; i < Q->nedge; i++) {
    o1 = Q->edge[2*i];
    o2 = Q->edge[2*i+1];
    n1 = Q->pop[o1];
    n2 = Q->pop[o2];
    
    c = gcd(n1, n2);
    
    d = random_factor (c);
    
    Q->weight[2*i] = n2/d;
    Q->weight[2*i+1] = n1/d;
    
  }

  return 0;
}

/**
 * @ingroup qgraph
 * @brief Use a complete quotient graph to create a full graph with the desired automorphism group
 * @param[in] Q The complete quotient graph
 * @param[out] G The graph to initialize
 * @return
 * @details
 *
 * @author Isaac Klickstein
 */
int SBD_qgraph_wire (SBD_qgraph_t *Q, igraph_t *G) {
  
  // Determine the number of nodes and edges to create in G
  int n = 0;
  int m = 0;
  
  // Sum the orbit populations
  for (int k = 0; k < Q->nnode; k++) {
    n += Q->pop[k];
  }
  
  // Add the intra-orbit edges
  for (int i = 0; i < Q->nnode; i++) {
    m += Q->loop[i] * Q->pop[i] / 2;
  }
  
  // Add the inter-orbit edges
  for (int i = 0; i < Q->nedge; i++) {
    m += Q->pop[Q->edge[2*i]] * Q->weight[2*i];
  }
  
  // Determine the first index of the nodes in each orbit
  int *first = (int*) calloc (Q->nnode, sizeof(int));
  first[0] = 0;
  for (int i = 1; i < Q->nnode; i++) {
    first[i] = first[i-1] + Q->pop[i-1];
  }
  
  
  // Initialize the edge weight
  igraph_vector_t edge;
  igraph_vector_init (&edge, 2*m);
  
  
  // Determine the intra-orbit edges
  int count = 0;
  int v1, v2;
  for (int k = 0; k < Q->nnode; k++) { // If there is an even number of edges...
    if (Q->loop[k] > 0) {
      if (Q->loop[k]%2 == 0) {
        for (int i = 0; i < Q->pop[k]; i++) {
          v1 = i;
          for (int j = 1; j <= Q->loop[k]/2; j++) {
            v2 = (i+j) % Q->pop[k];
            if (v1 < v2) {
              VECTOR(edge)[2*count]   = v1 + first[k];
              VECTOR(edge)[2*count+1] = v2 + first[k];
              count += 1;
            }
            v2 = (i+Q->pop[k]-j) % Q->pop[k];
            if (v1 < v2) {
              VECTOR(edge)[2*count]   = v1 + first[k];
              VECTOR(edge)[2*count+1] = v2 + first[k];
              count += 1;
            }
          }
        }
      }
      else { // otherwise if there is an odd number of edges
        for (int i = 0; i < Q->pop[k]; i++) {
          v1 = i;
          for (int j = 1; j <= (Q->loop[k]-1)/2; j++) {
            v2 = (i+j) % Q->pop[k];
            if (v1 < v2) {
              VECTOR(edge)[2*count] = v1 + first[k];
              VECTOR(edge)[2*count+1] = v2 + first[k];
              count += 1;
            }
            v2 = (i+Q->pop[k]-j) % Q->pop[k];
            if (v1 < v2) {
              VECTOR(edge)[2*count] = v1 + first[k];
              VECTOR(edge)[2*count+1] = v2 + first[k];
              count += 1;
            }
          }
          v2 = (i + Q->pop[k]/2) % Q->pop[k];
          if (v1 < v2) {
            VECTOR(edge)[2*count] = v1 + first[k];
            VECTOR(edge)[2*count+1] = v2 + first[k];
            count += 1;
          }
        }
      }
    }
  }
  
  int medge;
  int o1, o2, n1, n2;
  int f12, f21;
  int c, d1, d2;
  int a;
  
  // Determine the inter-orbit edges
  for (int k = 0; k < Q->nedge; k++) {
    o1 = Q->edge[2*k];   n1 = Q->pop[o1];
    o2 = Q->edge[2*k+1]; n2 = Q->pop[o2];
    f12 = Q->weight[2*k]; f21 = Q->weight[2*k+1];
    medge = n1*f12;
    c = gcd(n1, n2);
    d1 = n1/c;
    d2 = n2/c;
    
    a = c*f12/n2;
    
    /*
    printf("Summary for inter-orbit edges\n");
    printf("  Edge (%i,%i)\n", o1, o2);
    printf("  Pop(%i) = %i, Pop(%i) = %i\n", o1, n1, o2, n2);
    printf("  Weight(%i->%i) = %i\n", o2, o1, f12);
    printf("  Weight(%i->%i) = %i\n", o1, o2, f21);
    printf("  Num. of edges = %i\n", medge);
    printf("  c = GCD(%i,%i) = %i\n", n1, n2, c);
    printf("  d1 = %i, d2 = %i\n", d1, d2);
    printf("  Pattern length: %i\n", a);
    */
    
    int *b = (int*) calloc (a, sizeof(int));
    
    for (int i = 0; i < a; i++) {
      b[i] = 1;
    }
    for (int i = 0; i < c-a; i++) {
      b[rand() % a] += 1;
    }
    
    for (int i = 0; i < n1; i++) {
      v1 = i;
      for (int r1 = 0; r1 < d2; r1++) {
        for (int r2 = 1; r2 <= a; r2++) {
          v2 = i + r1*c;
          for (int j = 1; j <= r2; j++) {
            v2 += b[j-1];
          }
          v2 = v2 % n2;
          
          VECTOR(edge)[2*count] = v1 + first[o1];
          VECTOR(edge)[2*count+1] = v2 + first[o2];
          count += 1;
        }
      }
    }
    
    
    free (b);
  }
  
  // With the edges determined, create the graph
  igraph_create (G, &edge, n, IGRAPH_UNDIRECTED);

  // Free local memory
  free (first);
  igraph_vector_destroy (&edge);
  return 0;
}

