/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

//library to include
#include "DDalphaAMG.h"

enum { T, Z, Y, X };
#define printf0(...) do{if(rank==0)printf(__VA_ARGS__);}while(0)

int rank;
int iters=100;
int coarsest_local_lattice[4];
MPI_Comm comm_cart;
DDalphaAMG_init init;
DDalphaAMG_parameters params;
DDalphaAMG_status status;
char * options = "L:p:B:l:i:wmet:h:V:v";

/*
 * Setting standard values for DDalphaAMG_init
 */
void standard_init() {
  coarsest_local_lattice[T] = 1;
  coarsest_local_lattice[Z] = 1;
  coarsest_local_lattice[Y] = 1;
  coarsest_local_lattice[X] = 1;

  init.kappa = 0.142857143;
  init.mu = 0.;
  init.csw = 0.;

  init.bc = 1;

  init.number_of_levels = 3;

  init.procs[T] = 1;
  init.procs[Z] = 1;
  init.procs[Y] = 1;
  init.procs[X] = 1;

  init.number_openmp_threads = 1;

  init.block_lattice[T]=2;
  init.block_lattice[Z]=2;
  init.block_lattice[Y]=2;
  init.block_lattice[X]=2;

  init.comm_cart = MPI_COMM_WORLD;
  init.Cart_coords = NULL;
  init.Cart_rank = NULL;
  init.init_file = NULL;
}

/*
 * Printing implemented parameters
 */
void help( char * arg0 ) {
  static int printed = 0;
  if(!printed) {
    printf0("\n\n");
    printf0("Usage: %s [<option(s)>]\n", arg0);
    printf0("   -L T Z Y X   Local lattice size on the coarsest level in each direction\n");
    printf0("   -p T Z Y X   Processors in each direction\n");
    printf0("   -B T Z Y X   Block size in each direction on first level. A block size of (2 2 2 2) is used for the following levels.\n");
    printf0("   -l #         Number of levels, l (from 1 to 4)\n");
    printf0("   -i #         Number of iterations\n");
    printf0("   -w           use c_sw term (without the coarse self-coupling is diagonal)\n");
    printf0("   -m           use twisted mass term (without the coarse self-coupling is hermitian)\n");
    printf0("   -e           use ND twisted mass term (this double the size of the vectors)\n");
    printf0("   -t #         Number of OpenMp threads\n");
    printf0("   -V 1 [2] [3] Basis vectors between each level (l-1 numbers)\n");
    printf0("   -v           Verbose\n");
  }
  printf0("\n\n");
  printed++;
} 

/*
 * Printing implemented parameters
 */
void read_init_arg(int argc, char *argv[] ) {

  int opt, p=0, fail=0, mu;
  optind = 0;
  while ((opt = getopt(argc, argv, options)) != -1) {
    switch (opt) {
    case 'L':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -L.\n");
	  p++;
	  fail++;
	  break;
	}
	coarsest_local_lattice[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -L.\n");
	p++;
	fail++;
      }
      break;
    case 'p':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -p.\n");
	  p++;
	  fail++;
	  break;
	}
	init.procs[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -p.\n");
	p++;
      }
      break;
    case 'B':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -B.\n");
	  p++;
	  fail++;
	  break;
	}
	init.block_lattice[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -B.\n");
	p++;
      }
      break;
    case 'l':
      init.number_of_levels = atoi(optarg);
      break;
    case 'i':
      iters = atoi(optarg);
      break;
    case 'w':
      init.csw = 1.;
      break;
    case 'm':
      init.mu = 1.;
      break;
    case 't':
      init.number_openmp_threads = atoi(optarg);
      break;
    case '?':
    case 'h':
      p++;
      break;
    default: 
      break;
    }
  }
  
  if(p) {
    help(argv[0]);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  if(fail)
    MPI_Abort(MPI_COMM_WORLD,0);
}

void read_params_arg(int argc, char *argv[] ) {

  int opt, p=0, fail=0, mu;
  optind = 0;
  while ((opt = getopt(argc, argv, options)) != -1) {
    switch (opt) {
    case 'V':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 2) {
	  printf0("Error: too many arguments in -V.\n");
	  p++;
	  fail++;
	  break;
	}
	params.mg_basis_vectors[mu] = atoi(argv[optind]);
	mu++;        
      }
      break;
    case 'e':
      params.epsbar = 1;
      break;
    case 'v':
      params.print = 1;
      break;
    default: 
      break;
    }
  }

  if(p) {
    help(argv[0]);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  if(fail)
    MPI_Abort(MPI_COMM_WORLD,0);
}

int main( int argc, char *argv[] ) {
  
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  standard_init();
  read_init_arg(argc, argv);
  
  for(int mu=0; mu<4; mu++) {
    init.global_lattice[mu] = coarsest_local_lattice[mu]*init.procs[mu]*init.block_lattice[mu];
    for(int l=0; l<init.number_of_levels-2; l++) {
      init.global_lattice[mu] *= 2;
    }
  }
  printf0("global size: %d %d %d %d\n",init.global_lattice[T],init.global_lattice[X],init.global_lattice[Y],init.global_lattice[Z]);

  printf0("Running initialization...\n");
  DDalphaAMG_initialize( &init, &params, &status );
  printf0("Initialized %d levels in %.2f sec\n", status.success, status.time);

  int nlvl = status.success;
  read_params_arg(argc, argv);
  
  comm_cart =  DDalphaAMG_get_communicator();
  MPI_Comm_rank( comm_cart, &rank );

  // we don't want to waste time with the setup
  for(int l=0; l<nlvl; l++) {
    params.setup_iterations[l]=0;
    params.smoother_iterations=0;
  }
  
  printf0("Running updating\n");
  DDalphaAMG_update_parameters( &params, &status );
  if (status.success)
    printf0("Updating time %.2f sec\n", status.time);

  /*
   * Reading the configuration. In plaq, it returns the plaquette value
   *  if provided in the configuration header.
   */
  double *gauge_field;
  int vol = init.global_lattice[T] * init.global_lattice[X] * init.global_lattice[Y] *
    init.global_lattice[Z] / init.procs[T] / init.procs[X] / init.procs[Y] / init.procs[Z];
  gauge_field = (double *) calloc(18*4*vol,sizeof(double));

  printf0("Using unitary config config.\n");
  for(int i=0; i<4*vol; i++) {
    gauge_field[i*18+0] = 1.;
    gauge_field[i*18+8] = 1.;
    gauge_field[i*18+16] = 1.;
  }
      
  printf0("Setting config.\n");
  DDalphaAMG_set_configuration( gauge_field, &status );
  printf0("Setting configuration time %.2f sec\n", status.time);
  printf0("Computed plaquette %.13lf\n", status.info);
  free(gauge_field);

  printf0("Running setup\n");
  DDalphaAMG_setup( &status );
  printf0("Run %d setup iterations in %.2f sec (%.1f %% on coarse grid)\n", status.success,
	  status.time, 100.*(status.coarse_time/status.time));
  printf0("Total iterations on fine grid %d\n", status.iter_count);
  printf0("Total iterations on coarse grids %d\n", status.coarse_iter_count);

  /*
   * Defining fine and coarse vector randomly.
   */
  float *vector1[nlvl], *vector2[nlvl];
  
  for ( int i=0; i<nlvl; i++ ) {
    vector1[i] = DDalphaAMG_coarse_vector_alloc(i);
    vector2[i] = DDalphaAMG_coarse_vector_alloc(i);
    DDalphaAMG_coarse_vector_rand(vector1[i], i);
  }

  for ( int i=1; i<nlvl; i++ ) {
    printf0("\nTesting RP=1 on level %d\n",i);
    DDalphaAMG_prolongate(vector2[i-1], vector1[i], i-1, &status);
    DDalphaAMG_restrict(vector2[i], vector2[i-1], i-1, &status);

    float res = DDalphaAMG_coarse_vector_residual(vector2[i], vector1[i], i);
    printf0("Result (1-RP)v = %e\n", res);
  }

  for ( int i=1; i<nlvl; i++ ) {
    printf0("\nTesting coarse operator on level %d\n",i);
    DDalphaAMG_apply_coarse_operator(vector2[i], vector1[i], i, &status);

    DDalphaAMG_prolongate(vector1[i-1], vector1[i], i-1, &status);
    DDalphaAMG_apply_coarse_operator(vector2[i-1], vector1[i-1], i-1, &status);
    DDalphaAMG_restrict(vector1[i], vector2[i-1], i-1, &status);

    float res = DDalphaAMG_coarse_vector_residual(vector2[i], vector1[i], i);
    printf0("Result (D_c-RDP)v = %e\n", res);

    // restoring vector1[i]
    DDalphaAMG_restrict(vector1[i], vector1[i-1], i-1, &status);
  }

  printf0("\nStarting bechmarks:\n");
  for ( int i=1; i<nlvl; i++ ) {
    double time = 0;
    for ( int j=0; j<iters; j++ ) {
      DDalphaAMG_prolongate(vector1[i-1], vector1[i], i-1, &status);
      time += status.time;
    }
    printf0("Level %d, prolongation from level %d to level %d: %e averaged over %d iters\n", i, i, i-1, time/iters, iters);
    time = 0;
    for ( int j=0; j<iters; j++ ) {
      DDalphaAMG_restrict(vector1[i], vector1[i-1], i-1, &status);
      time += status.time;
    }
    printf0("Level %d, restriction from level %d to level %d: %e averaged over %d iters\n", i, i-1, i, time/iters, iters);
    time = 0;
    for ( int j=0; j<iters; j++ ) {
      DDalphaAMG_apply_coarse_operator(vector1[i], vector1[i], i, &status);
      time += status.time;
    }
    printf0("Level %d, coarse operator on level %d: %e averaged over %d iters\n", i, i, time/iters, iters);
  }

  
  for ( int i=0; i<nlvl; i++ ) {
    DDalphaAMG_coarse_vector_free(vector1[i], i);
    DDalphaAMG_coarse_vector_free(vector2[i], i);
  }
  
  DDalphaAMG_finalize();
  MPI_Finalize();
  return 0;
}
