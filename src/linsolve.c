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

#include "main.h"
#include "linsolve.h"

void fgmres_MP_struct_init( gmres_MP_struct *p ) {
  fgmres_float_struct_init( &(p->sp) );
  fgmres_double_struct_init( &(p->dp) );
}


void fgmres_MP_struct_alloc( int m, int n, const int vl_type, double tol, const int prec_kind, 
                             void (*precond)(), gmres_MP_struct *p, level_struct *l ) {
  long int total=0; 
  int i, k=0, n_vl=g.num_rhs_vect;
  
  p->dp.restart_length = m;                      p->sp.restart_length = m;           
  p->dp.num_restart = n;                         p->sp.num_restart = n;
  p->dp.preconditioner = NULL;                   p->sp.preconditioner = precond;
  if ( g.method == 6 ) {
  p->dp.eval_operator = g5D_plus_clover_double;  p->sp.eval_operator = g5D_plus_clover_float;
  } else {
  p->dp.eval_operator = d_plus_clover_double_new; p->sp.eval_operator = d_plus_clover_float_new;
  }
  p->dp.tol = tol;                               p->sp.tol = MAX(tol,1E-5);
  p->dp.kind = _NOTHING;                         p->sp.kind = prec_kind;
  p->dp.timing = 1;                              p->sp.timing = 1;
                                 
  p->dp.print = g.vt.evaluation?0:1;             p->sp.print = g.vt.evaluation?0:1;
  p->dp.initial_guess_zero = 1;                  p->sp.initial_guess_zero = 1;
  p->dp.v_start = 0;                             p->sp.v_start = 0;
  p->dp.v_end = l->inner_vector_size;            p->sp.v_end = l->inner_vector_size;
  
  p->dp.op = &(g.op_double);                     p->sp.op = &(l->s_float.op);
  
  g.p.op = &(g.op_double);
  if ( g.method == 6 ) {
    g.p.eval_operator = g5D_plus_clover_double;
  } else {
    g.p.eval_operator = d_plus_clover_double_new;
  }
  
#ifdef HAVE_TM1p1
  n_vl*=2;
#endif

  // double precision part
  total = 0;
  total += (m+1)*m*n_vl; // Hessenberg matrix
  MALLOC( p->dp.H, complex_double*, m );
  total += 4*(m+1)*n_vl; // y, gamma, c, s
  p->dp.total_storage = total;
  // precomputed storage amount
  
  p->dp.H[0] = NULL; // allocate connected memory
  MALLOC( p->dp.H[0], complex_double, total );
  
  // reserve storage
  total = 0;
  // H
  for ( i=1; i<m; i++ )
    p->dp.H[i] = p->dp.H[0] + i*(m+1)*n_vl;
  total += m*(m+1)*n_vl;
  // y
  p->dp.y = p->dp.H[0] + total; total += (m+1)*n_vl;
  // gamma
  p->dp.gamma = p->dp.H[0] + total; total += (m+1)*n_vl;
  // c
  p->dp.c = p->dp.H[0] + total; total += (m+1)*n_vl;
  // s
  p->dp.s = p->dp.H[0] + total; total += (m+1)*n_vl;
  // x
  vector_double_alloc( &(p->dp.x), vl_type, n_vl, l, no_threading );
  // r
  vector_double_alloc( &(p->dp.r), vl_type, n_vl, l, no_threading );
  // b
  vector_double_alloc( &(p->dp.b), vl_type, n_vl, l, no_threading );
  
  ASSERT( p->dp.total_storage == total );
  
  
  // single precision part
  total = 0;
  MALLOC( p->sp.V, vector_float, m+1 );
  if ( precond != NULL ) {
    if ( prec_kind == _RIGHT ) {
      k = m+1;
    } else {
      k = 1;
    }
    MALLOC( p->sp.Z, vector_float, k );
  }
  p->sp.total_storage = total;
  // precomputed storage amount
  
  // reserve storage
  total = 0;
  // w
  vector_float_alloc( &(p->sp.w), vl_type, n_vl, l, no_threading );
  // V 
  for ( i=0; i<m+1; i++ ) {
    vector_float_init( &(p->sp.V[i]) );
    vector_float_alloc( &(p->sp.V[i]), vl_type, n_vl, l, no_threading );
  }
  // Z
  if ( precond != NULL ) {
    for ( i=0; i<k; i++ ) {
      vector_float_init( &(p->sp.Z[i]) );
      vector_float_alloc( &(p->sp.Z[i]), vl_type, n_vl, l, no_threading );
    }
  }
  
  ASSERT( p->sp.total_storage == total );
}  
   
   
void fgmres_MP_struct_free( gmres_MP_struct *p, level_struct *l ) {
   
  // single precision
  vector_float_free( &(p->sp.w), l, no_threading );
  FREE( p->sp.V, vector_float, p->sp.restart_length+1 );
  if ( p->sp.Z != NULL )
    FREE( p->sp.Z, vector_float, p->sp.kind==_RIGHT?p->sp.restart_length+1:1 );
  
  // double precision
  FREE( p->dp.H[0], complex_double, p->dp.total_storage );
  FREE( p->dp.H, complex_double*, p->dp.restart_length );
  vector_double_free( &(p->dp.x), l, no_threading );
  vector_double_free( &(p->dp.r), l, no_threading );
  vector_double_free( &(p->dp.b), l, no_threading );
} 
  
  
int fgmres_MP( gmres_MP_struct *p, level_struct *l, struct Thread *threading ) {
  
/*********************************************************************************
* Uses FGMRES to solve the system D x = b, where b is taken from p->b and x is 
* stored in p->x.                                                              
*********************************************************************************/  
  
  // start and end indices for vector functions depending on thread
  ASSERT( g.mixed_precision );
  
  int start;
  int end;
  
  int j=-1, finish=0, iter=0, il, ol, n_vect=g.num_rhs_vect, i, jj;
  complex_double gamma0[n_vect];//gamma0=0;
  double beta[n_vect]; //beta=0;

  double t0=0, t1=0;
  double norm_r0[n_vect], gamma_jp1[n_vect], gamma0_real[n_vect], gamma_tot, H_tot, gamma_tot2;//norm_r0=1, gamma_jp1=1
  complex_float gamma_float[n_vect];
  
  VECTOR_LOOP(i, n_vect, jj, norm_r0[i+jj]=1;
                             gamma_jp1[i+jj]=1;)
  
  START_LOCKED_MASTER(threading)
#ifndef WILSON_BENCHMARK
  if ( l->depth==0 && ( p->dp.timing || p->dp.print ) ) prof_init( l );
#endif
  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->dp.tol = g.coarse_tol;
  if ( l->depth > 0 ) p->dp.timing = 1;
  if ( l->depth == 0 ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( p->dp.print && g.print > 0 ) printf0("+----------------------------------------------------------+\n");
#endif
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  //compute_core_start_end(p->dp.v_start, p->dp.v_end, &start, &end, l, threading);
  
  // Outer loop in double precision
  for( ol=0; ol<p->dp.num_restart && finish==0; ol++ )  {
    
    if( ol == 0 && p->dp.initial_guess_zero ) {
      //vector_double_copy( &(p->dp.r), &(p->dp.b), start, end, l );
      vector_double_copy_new( &(p->dp.r), &(p->dp.b), l, threading );
    } else {
      apply_operator_double( &(p->dp.r), &(p->dp.x), &(p->dp), l, threading ); // compute r <- D*x
      //vector_double_minus( &(p->dp.r), &(p->dp.b), &(p->dp.r), start, end, l ); // compute r <- b - r
      vector_double_minus_new( &(p->dp.r), &(p->dp.b), &(p->dp.r), l, threading );
    }
    //gamma0 = (complex_double) global_norm_double( &(p->dp.r), p->dp.v_start, p->dp.v_end, l, threading ); // gamma_0 = norm(r)
    global_norm_double_new( gamma0_real, &(p->dp.r), l, threading );
    VECTOR_LOOP(i, n_vect, jj, gamma0[i+jj]=gamma0_real[i+jj];)

    START_MASTER(threading)
    //p->dp.gamma[0] = gamma0;
    VECTOR_LOOP(i, n_vect, jj, p->dp.gamma[i+jj] = gamma0[i+jj];)
    
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if( ol == 0) {
     if (l->depth == 0 && !p->dp.initial_guess_zero) {
       //norm_r0 = global_norm_double( &(p->dp.b), start, end, l, threading );
       global_norm_double_new( norm_r0, &(p->dp.b), l, threading );
       for( i=0; i<n_vect; i++ ) 
         printf0("| initial guess relative residual (%d):            %le |\n", i, creal(gamma0[i])/norm_r0[i]);
     } else {
       //norm_r0 = creal(gamma0);
       VECTOR_LOOP(i, n_vect, jj, norm_r0[i+jj]= creal(gamma0[i+jj]);)
     }
    }
/*#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    else {
      if ( p->dp.print && g.print > 0 ) {
        START_MASTER(threading)
        printf0("+----------------------------------------------------------+\n");
        for( i=0; i<n_vect; i++ ) 
          printf0("| restarting ...          true residual norm (%d): %6e |\n", i, creal(gamma0[i])/norm_r0[i] );
        printf0("+----------------------------------------------------------+\n");
        END_MASTER(threading)
      }
    }
#endif*/
    trans_float_new( &(p->sp.V[0]), &(p->dp.r), l->s_float.op.translation_table, l, threading );
    //vector_float_real_scale( &(p->sp.V[0]), &(p->sp.V[0]), (float)(1/p->dp.gamma[0]), start, end, l ); // V[0] <- r / gamma_0
    VECTOR_LOOP(i, n_vect, jj, gamma_float[i+jj]= (complex_float) p->dp.gamma[0*n_vect+i+jj];)
    vector_float_real_scale_new( &(p->sp.V[0]), &(p->sp.V[0]), gamma_float, 0, 1, l, threading );
    // inner loop in single precision
    for( il=0; il<p->dp.restart_length && finish==0; il++) {
      j = il; iter++;
      arnoldi_step_MP_new( p->sp.V, p->sp.Z, &(p->sp.w), p->dp.H, p->dp.y, j, p->sp.preconditioner, &(p->sp), l, threading );
      H_tot=0;
      VECTOR_LOOP(i, n_vect, jj, H_tot += cabs( p->dp.H[j][(j+1)*n_vect+i+jj] );)
      //if ( cabs( p->dp.H[j][j+1] ) > 1E-15 )
      if ( H_tot > n_vect*1E-15 ) {
        qr_update_double( p->dp.H, p->dp.s, p->dp.c, p->dp.gamma, j, l, threading );
        //gamma_jp1 = cabs( p->dp.gamma[j+1] );          
        VECTOR_LOOP(i, n_vect, jj, gamma_jp1[i+jj] = cabs( p->dp.gamma[(j+1)*n_vect+i+jj] );)

        if ( iter%10 == 0 || p->sp.preconditioner != NULL || l->depth > 0 ) {
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
          START_MASTER(threading)
          if ( p->sp.print && g.print > 0 )
            for( i=0; i<n_vect; i++ )
              printf0("| vector %d, approx. rel. res. after  %-6d iterations: %e |\n", i, iter, gamma_jp1[i]/norm_r0[i] );
          END_MASTER(threading)
#endif
        }
        gamma_tot=0;
        VECTOR_LOOP(i, n_vect, jj, gamma_tot += gamma_jp1[i+jj]/norm_r0[i+jj];)

        //if( gamma_jp1/norm_r0 < p->dp.tol || gamma_jp1/norm_r0 > 1E+5 )  // if satisfied ... stop
        if( gamma_tot < n_vect*p->dp.tol || gamma_tot > n_vect*1E+5 ) {
          finish = 1;
          START_MASTER(threading)
            if ( gamma_tot > n_vect*1E+5 ) printf0("Divergence of fgmres_MP, iter = %d, level=%d\n", iter, l->level );
          END_MASTER(threading)
        }
        gamma_tot2=0;
        VECTOR_LOOP(i, n_vect, jj, gamma_tot2 += gamma_jp1[i+jj]/creal(gamma0[i+jj]);)
        //if( gamma_jp1/creal(gamma0) < p->sp.tol )
        if( gamma_tot2 < n_vect*p->sp.tol ){  
          break;
        }
      } else {
        finish = 1;
      }
    } // end of a single restart
    compute_solution_MP_new( &(p->sp.w), (p->sp.preconditioner&&p->sp.kind==_RIGHT)?p->sp.Z:p->sp.V,
                         p->dp.y, p->dp.gamma, p->dp.H, j, &(p->sp), l, threading );
                           
    trans_back_float_new( &(p->dp.r), &(p->sp.w), l->s_float.op.translation_table, l, threading );
    if ( ol == 0 ) {
      //vector_double_copy( &(p->dp.x), &(p->dp.r), start, end, l );
      vector_double_copy_new(&(p->dp.x), &(p->dp.r), l, threading);
    } else {
      //vector_double_plus( &(p->dp.x), &(p->dp.x), &(p->dp.r), start, end, l );
      vector_double_plus_new( &(p->dp.x), &(p->dp.x), &(p->dp.r), l, threading );
    }
  } // end of fgmres
  
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.total_time = t1-t0; g.iter_count = iter; g.norm_res = gamma_tot; }
  END_LOCKED_MASTER(threading)
  
  if ( p->dp.print ) {
#ifdef FGMRES_RESTEST
    apply_operator_double( &(p->dp.r), &(p->dp.x), &(p->dp), l, threading );
    //vector_double_minus( &(p->dp.r), &(p->dp.b), &(p->dp.r), start, end, l );
    vector_double_minus_new( &(p->dp.r), &(p->dp.b), &(p->dp.r), l, threading ); 
    //beta = global_norm_double( &(p->dp.r), p->dp.v_start, p->dp.v_end, l, threading );
    global_norm_double_new( beta, &(p->dp.r), l, threading );
#else
    VECTOR_LOOP(i, n_vect, jj, beta[i+jj] = creal(gamma_jp1[i+jj]);)
#endif
    START_MASTER(threading)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if ( g.print > 0 ) printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|    FGMRES MP iterations: %-6d coarse average: %-6.2lf   |\n", iter,
            ((double)g.coarse_iter_count)/((double)iter) );
    for( i=0; i<n_vect; i++ )
      printf0("| exact relative residual %d: ||r||/||b|| = %e    |\n",i, beta[i]/norm_r0[i] );
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", t1-t0 );
    if ( g.coarse_time > 0 ) 
      printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.coarse_time, 100*(g.coarse_time/(t1-t0)) );
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }

  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", iter );
      printf0("solve time: %le seconds\n", t1-t0 );
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
    END_LOCKED_MASTER(threading)
    }
  }

  if ( l->depth == 0 && ( p->dp.timing || p->dp.print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
#ifndef WILSON_BENCHMARK
    prof_print( l );
#endif
    END_MASTER(threading)
  }
  
  return iter;
}


void arnoldi_step_MP( vector_float *V, vector_float *Z, vector_float *w,
                      complex_double **H, complex_double* buffer, int j, void (*prec)(),
                      gmres_float_struct *p, level_struct *l, struct Thread *threading ) {
  
  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
  int i;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  
  if ( prec != NULL ) {
    if ( p->kind == _LEFT ) {
      apply_operator_float( &Z[0], &V[j], p, l, threading );
      prec( w, NULL, &Z[0], _NO_RES, l, threading );
    } else {
      if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 ) ) {
        prec( &Z[j], w, &V[j], _NO_RES, l, threading );
        // obtains w = D * Z[j] from Schwarz
      } else {
        prec( &Z[j], NULL, &V[j], _NO_RES, l, threading );
        apply_operator_float( w, &Z[j], p, l, threading ); // w = D*Z[j]
      }
    }
  } else {
    apply_operator_float( w, &V[j], p, l, threading ); // w = D*V[j]
  }

  complex_double tmp[j+1];
  process_multi_inner_product_MP( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    buffer[i] = tmp[i];
  
  if ( g.num_processes > 1 ) {
    PROF_double_START( _ALLR );
    MPI_Allreduce( buffer, H[j], j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_float.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  } else {
    for( i=0; i<=j; i++ )
      H[j][i] = buffer[i];
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // orthogonalization
  complex_float alpha[j+1];
  for( i=0; i<=j; i++ )
    alpha[i] = (complex_float) -H[j][i];
  vector_float_multi_saxpy( w, V, alpha, 1, j+1, start, end, l );
  
  complex_double tmp2 = global_norm_MP( w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  H[j][j+1] = tmp2;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // V_j+1 = w / H_j+1,j
  if ( cabs_double( H[j][j+1] ) > 1e-15 )
    vector_float_real_scale( &V[j+1], w, (float)(1/H[j][j+1]), start, end, l );
}


void arnoldi_step_MP_new( vector_float *V, vector_float *Z, vector_float *w,
                      complex_double **H, complex_double* buffer, int j, void (*prec)(),
                      gmres_float_struct *p, level_struct *l, struct Thread *threading ) {
  
  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
  int i, n_vect=g.num_rhs_vect, n, jj;
  double H_tot;
  complex_float H_float[n_vect];
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  //compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  
  if ( prec != NULL ) {
    if ( p->kind == _LEFT ) {
      apply_operator_float( &Z[0], &V[j], p, l, threading );
      prec( w, NULL, &Z[0], _NO_RES, l, threading );
    } else {
      if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 ) ) {
        prec( &Z[j], w, &V[j], _NO_RES, l, threading );
        // obtains w = D * Z[j] from Schwarz
      } else {
        prec( &Z[j], NULL, &V[j], _NO_RES, l, threading );
        apply_operator_float( w, &Z[j], p, l, threading ); // w = D*Z[j]
      }
    }
  } else {
    apply_operator_float( w, &V[j], p, l, threading ); // w = D*V[j]
  }

  complex_double tmp[(j+1)*n_vect];
  process_multi_inner_product_MP_new( j+1, tmp, V, w, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(n, n_vect, jj, buffer[i*n_vect+n+jj] = tmp[i*n_vect+n+jj];)

  if ( g.num_processes > 1 ) {
    PROF_double_START( _ALLR );
    MPI_Allreduce( buffer, H[j], (j+1)*n_vect, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_float.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  } else {
    for( i=0; i<=j; i++ )
      VECTOR_LOOP(n, n_vect, jj, H[j][i*n_vect+n+jj] = buffer[i*n_vect+n+jj];)
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  complex_float alpha[(j+1)*n_vect]; 
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(n, n_vect, jj, alpha[i*n_vect+n+jj] = (complex_float) H[j][i*n_vect+n+jj];)
  for( i=0; i<=j; i++ )
    vector_float_saxpy_new( w, w, &V[i], alpha, i, -1, l, threading );
  /*// orthogonalization
  complex_float alpha[(j+1)*n_vect];
 
  for( i=0; i<=j; i++ )
    for( n_vec=0; n_vec<n_vect; n_vec++ )
      alpha[i*n_vect+n_vec] = (complex_float) -H[j][i*n_vect+n_vec];
  vector_float_multi_saxpy_new( w, V, alpha, 1, j+1, l, threading );
  */
  double tmp2[n_vect];
  global_norm_MP_new( tmp2, w, l, threading );
  START_MASTER(threading)
  VECTOR_LOOP(n, n_vect, jj, H[j][(j+1)*n_vect+n+jj] = tmp2[n+jj];)

  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // V_j+1 = w / H_j+1,j
  H_tot=0;
  VECTOR_LOOP(n, n_vect, jj, H_tot += cabs_double( H[j][(j+1)*n_vect+n+jj] );)
  
  if ( H_tot > n_vect*1e-15 ){
    VECTOR_LOOP(n, n_vect, jj, H_float[n+jj]= (complex_float) H[j][(j+1)*n_vect+n+jj];)
   vector_float_real_scale_new( &V[j+1], w, H_float, 0, 1, l, threading );
  }
}


void compute_solution_MP( vector_float *x, vector_float *V, complex_double *y,
                          complex_double *gamma, complex_double **H, int j,
                          gmres_float_struct *p, level_struct *l, struct Thread *threading ) {
  
  int i, k;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  START_MASTER(threading)
  
  PROF_double_START( _SMALL2 );
  
  // backward substitution
  for ( i=j; i>=0; i-- ) {
    y[i] = gamma[i];
    for ( k=i+1; k<=j; k++ ) {
      y[i] -= H[k][i]*y[k];
    }
    y[i] /= H[i][i];
  }
  
  PROF_double_STOP( _SMALL2, ((j+1)*(j+2))/2 + j+1 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // x = V*y
  vector_float_scale( x, &V[0], (complex_float) y[0], start, end, l );

  complex_float alpha[j];
  for ( i=1; i<=j; i++ )
    alpha[i-1] = (complex_float) y[i];
  vector_float_multi_saxpy( x, &V[1], alpha, 1, j, start, end, l );
}


void compute_solution_MP_new( vector_float *x, vector_float *V, complex_double *y,
                          complex_double *gamma, complex_double **H, int j,
                          gmres_float_struct *p, level_struct *l, struct Thread *threading ) {
  
  int i, k, n, jj, n_vect=g.num_rhs_vect;
  complex_float y_float[n_vect];
  // start and end indices for vector functions depending on thread
  //int start;
  //int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  //compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  START_MASTER(threading)
  
  PROF_double_START( _SMALL2 );
  
  // backward substitution
  for ( i=j; i>=0; i-- ) {
    VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] = gamma[i*n_vect+n+jj];)
    for ( k=i+1; k<=j; k++ ) {
      for ( n=0; n<n_vect; n++ )
        y[i*n_vect+n] -= H[k][i*n_vect+n]*y[k*n_vect+n];
    }
    VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] /= H[i][i*n_vect+n+jj];)
  }
  
  PROF_double_STOP( _SMALL2, ((j+1)*(j+2))/2 + j+1 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // x = V*y
  VECTOR_LOOP(n, n_vect, jj, y_float[n+jj]= (complex_float) y[0*n_vect+n+jj];)
  vector_float_scale_new( x, &V[0], y_float, 0, l, threading );

  complex_float alpha[j*n_vect];
  for ( i=1; i<=j; i++ )
    VECTOR_LOOP(n, n_vect, jj, alpha[i*n_vect+n+jj] = (complex_float) y[i*n_vect+n+jj];)
  for ( i=1; i<=j; i++ )
    vector_float_saxpy_new( x, x, &V[i], alpha, i, 1, l, threading );
//vector_float_multi_saxpy_new( x, &V[1], alpha, 1, j, l, threading );

}

