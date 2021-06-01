
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "SBD.h"

int main (int argc, char *argv[]) {
  
  // PARSE INPUTS
  bool laplacian = 0;
  bool use_orbits = 1;
  bool write_P = 0;
  bool write_B = 0;
  
  int opt = 1;
  
  while ((opt = getopt(argc, argv, "laPBOM")) != -1) {
    switch (opt) {
      case 'l':
        laplacian = 1;
        break;
      case 'a':
        laplacian = 0;
        break;
      case 'O':
        use_orbits = 1;
        break;
      case 'M':
        use_orbits = 0;
        break;
      case 'P':
        write_P = 1;
        break;
      case 'B':
        write_B = 1;
        break;
      default: 
        fprintf (stderr, "USAGE: ./runSBD <options> <directory>");
        exit (EXIT_FAILURE);
    }
  }
  
  if (optind >= argc) {
    fprintf (stderr, "Expected directory with graph supplied\n");
    exit (EXIT_FAILURE);
  }
  
  
  SBD_t SBD;
  igraph_t graph;
  char directory[128];
  char filename[256];

  sprintf (directory, "%s", argv[optind]);
  sprintf (filename, "%s/graph", directory);
  
  FILE *fin = fopen(filename, "r");
  
  if (fin == NULL) {
    fprintf (stderr, "The graph file was not supplied in the directory");
    exit (EXIT_FAILURE);
  }
  
  igraph_read_graph_edgelist (&graph, fin, 0, IGRAPH_UNDIRECTED);
  fclose (fin);

  // CREATE THE SBD INSTANCE
  SBD_from_graph (&SBD, &graph, use_orbits, laplacian);
  SBD_print (&SBD);
  
  int n = igraph_vcount (&graph);
  int m = igraph_ecount (&graph);
  
  // COMPUTE THE ORBITS 
  igraph_vector_int_t orbits;
  igraph_vector_int_init (&orbits, n);
  SBD_nauty_solve (&graph, &orbits);
  printf("Orbits: \n");
  igraph_vector_int_print (&orbits);
  
  // WRITE THE ORBITS TO A FILE
  sprintf (filename, "%s/orbits", directory);
  FILE *fout = fopen(filename, "w");

  for (int k = 0; k < n; k++) {
    fprintf(fout, "%i\n", VECTOR(orbits)[k]);
  }
  
  fclose (fout);
  
  // WRITE THE MBC TO A FILE
  sprintf (filename, "%s/MBC", directory);
  fout = fopen(filename, "w");
  SBD_belykh (&graph, &orbits);
  
  printf("MBC:\n");
  igraph_vector_int_print (&orbits);
  
  for (int k = 0; k < n; k++) {
    fprintf (fout, "%i\n", VECTOR(orbits)[k]);
  }
  
  fclose (fout);
  
  
  igraph_destroy (&graph);

  
  // COMPUTE THE BLOCK DIAGONALIZING TRANSFORMATION P
  igraph_matrix_t P;
  igraph_matrix_init (&P, n, n);
  
  double epsilon = 1E-13;
  
  SBD_lanczos (&SBD, epsilon, &P);
  
  if (write_P) {
  
    sprintf (filename, "%s/P", directory);
    fout = fopen (filename, "w");
  
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        fprintf (fout, "%E ", MATRIX(P,i,j));
      }
      fprintf (fout, "\n");
    }
    fclose (fout);
  }
  
  if (write_B) {  
    for (int k = 0; k < SBD.M; k++) {
      sprintf (filename, "%s/B_%i", directory, k);
      fout = fopen (filename, "w");
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          fprintf (fout, "%E ", MATRIX(SBD.B[k],i,j));
        }
        fprintf (fout, "\n");
      }
      fclose (fout);
    }
  }
  
  igraph_vector_int_t blocks;
  igraph_vector_int_init (&blocks, n);
  SBD_block_structure (&SBD, &blocks);
  
  igraph_vector_int_print (&blocks);
  
  SBD_free (&SBD);
  igraph_matrix_destroy (&P);
  igraph_vector_int_destroy (&orbits);
  igraph_vector_int_destroy (&blocks);
  
  
  return 0;
}
