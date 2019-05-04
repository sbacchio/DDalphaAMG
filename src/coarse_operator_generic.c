/*
 * copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
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

void coarse_operator_PRECISION_alloc( level_struct *l ) {
  
  int nd = l->next_level->num_inner_lattice_sites,
      k = l->next_level->num_parent_eig_vect*2;  
  l->next_level->D_size = k*k*4*nd;
  l->next_level->clover_size = ((k*(k+1))/2)*nd;
  l->next_level->block_size = ((k/2*(k/2+1)))*nd;
  
  operator_PRECISION_alloc( &(l->next_level->op_PRECISION), _ORDINARY, l->next_level );

}

void coarse_operator_PRECISION_free( level_struct *l ) {
  
  operator_PRECISION_free( &(l->next_level->op_PRECISION), _ORDINARY, l->next_level );
  
}

void coarse_operator_PRECISION_setup( vector_PRECISION *V, level_struct *l ) {
  
  double t0, t1;
  t0 = MPI_Wtime();
  
  vector_PRECISION buffer1, buffer2;
  buffer1.vector_buffer = l->vbuf_PRECISION[4].vector_buffer; buffer2.vector_buffer = l->vbuf_PRECISION[5].vector_buffer;
  
  int mu, n = l->num_eig_vect, i, j,
    D_size = l->next_level->D_size,
    clover_size = l->next_level->clover_size,
    block_size = l->next_level->block_size;
  void (*aggregate_self_coupling)() = (l->depth==0)?d_plus_clover_aggregate_PRECISION:coarse_aggregate_self_couplings_PRECISION,
       (*aggregate_neighbor_coupling)() = (l->depth==0)?d_neighbor_aggregate_PRECISION:coarse_aggregate_neighbor_couplings_PRECISION;
  void (*aggregate_block)() = (l->depth==0)?diagonal_aggregate_PRECISION:coarse_aggregate_block_diagonal_PRECISION;
   
  operator_PRECISION_define( &(l->next_level->op_PRECISION), l->next_level );
    
  for ( j=0; j<D_size; j++ )
    l->next_level->op_PRECISION.D[j] = _COMPLEX_PRECISION_ZERO;
  for ( j=0; j<clover_size; j++ )
    l->next_level->op_PRECISION.clover[j] = _COMPLEX_PRECISION_ZERO;
  for ( j=0; j<block_size; j++ )
    l->next_level->op_PRECISION.odd_proj[j] = _COMPLEX_PRECISION_ZERO;
  
  // for all test vectors V[i]:
  for ( i=0; i<n; i++ ) {
    for ( mu=0; mu<4; mu++ ) {
      // update ghost cells of V[i]
      negative_sendrecv_PRECISION( &V[i], mu, &(l->s_PRECISION.op.c), l );
    }
    // apply self coupling of block-and-2spin-restricted dirac operator for each aggregate
    aggregate_self_coupling( &buffer1, &buffer2, &V[i], &(l->s_PRECISION), l );
    // calculate selfcoupling entries of the coarse grid operator
    set_coarse_self_coupling_PRECISION( &buffer1, &buffer2, V, i, l );
    //odd_proj
    aggregate_block( &buffer1, &buffer2, &V[i], l->s_PRECISION.op.odd_proj, l );
    set_block_diagonal_PRECISION( &buffer1, &buffer2, V, i, l->next_level->op_PRECISION.odd_proj, l );
 
    for ( mu=0; mu<4; mu++ ) {
      // finish updating ghostcells of V[i]
      negative_wait_PRECISION( mu, &(l->s_PRECISION.op.c), l );      
      // apply 2spin-restricted dirac operator for direction mu for all aggregates
      aggregate_neighbor_coupling( &buffer1, &buffer2, &V[i], mu, &(l->s_PRECISION), l );      
      set_coarse_neighbor_coupling_PRECISION( &buffer1, &buffer2, V, mu, i, l );
    }
  }
  
  coarse_operator_PRECISION_setup_finalize( l, no_threading );

  t1 = MPI_Wtime();
  if ( g.print > 0 ) printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds\n", l->depth, t1-t0 );

}

void coarse_operator_PRECISION_setup_finalize( level_struct *l, struct Thread *threading ) {

  int block_size = l->next_level->block_size;
  
  l->next_level->op_PRECISION.m0 = l->s_PRECISION.op.m0;
#ifdef HAVE_TM    
  //tm_term
  PRECISION mf = (g.mu_factor[l->depth]) ? g.mu_factor[l->next_level->depth]/g.mu_factor[l->depth]:0;
  if ( mf*l->s_PRECISION.op.mu + mf*l->s_PRECISION.op.mu_even_shift == 0 &&
       mf*l->s_PRECISION.op.mu + mf*l->s_PRECISION.op.mu_odd_shift == 0 )
    buffer_PRECISION_define( l->next_level->op_PRECISION.tm_term, _COMPLEX_double_ZERO, 0, block_size, l->next_level );
  else
    tm_term_PRECISION_setup( mf*l->s_PRECISION.op.mu, mf*l->s_PRECISION.op.mu_even_shift,
                             mf*l->s_PRECISION.op.mu_odd_shift, &(l->next_level->op_PRECISION),
                             l->next_level, threading ); 
#endif
#ifdef HAVE_TM1p1
  //eps_term
  PRECISION ef = (g.epsbar_factor[l->depth]) ? g.epsbar_factor[l->next_level->depth]/g.epsbar_factor[l->depth]:0; 
  if ( ef*l->s_PRECISION.op.epsbar == 0 &&  ef*l->s_PRECISION.op.epsbar_ig5_even_shift == 0 &&
       ef*l->s_PRECISION.op.epsbar_ig5_odd_shift == 0 )
    buffer_PRECISION_define( l->next_level->op_PRECISION.epsbar_term, _COMPLEX_double_ZERO, 0, block_size, l->next_level );
  else
    epsbar_term_PRECISION_setup( ef*l->s_PRECISION.op.epsbar, ef*l->s_PRECISION.op.epsbar_ig5_even_shift,
                                 ef*l->s_PRECISION.op.epsbar_ig5_odd_shift, &(l->next_level->op_PRECISION),
                                 l->next_level, threading );
#endif

}

void set_block_diagonal_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, 
                                   vector_PRECISION *V, const int n, config_PRECISION block, level_struct *l ) {
  
  // U(x) = [ A 0      , A=A*, D=D*
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, columnwise
  // suitable for tm_term and odd_proj
  
  int i, j, k, m, k1, k2, num_aggregates = l->is_PRECISION.num_agg,
    num_eig_vect = l->next_level->num_parent_eig_vect,
    aggregate_size = l->num_inner_lattice_sites*l->num_parent_eig_vect*2/num_aggregates, 
    offset = l->num_parent_eig_vect,
    block_site_size = (num_eig_vect*(num_eig_vect+1));
  buffer_PRECISION spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_PRECISION block_pt;

  for ( k=0; k<=n; k++ ) {
    k1 = (n*(n+1))/2+k; k2 = (n*(n+1))/2+k+block_site_size/2;
  
    for ( j=0; j<num_aggregates; j++ ) {
      spin_0_1_pt = spin_0_1->vector_buffer + j*aggregate_size;
      spin_2_3_pt = spin_2_3->vector_buffer + j*aggregate_size;
      interpolation_data = V[k].vector_buffer + j*aggregate_size;
      block_pt = block + j*block_site_size;
      
      for ( i=0; i<aggregate_size; ) {
        // A
        for ( m=0; m<offset; m++, i++ )
          block_pt[ k1 ] += conj_PRECISION( interpolation_data[i] ) * spin_0_1_pt[i];
        // D
        for ( m=0; m<offset; m++, i++ )
          block_pt[ k2 ] += conj_PRECISION( interpolation_data[i] ) * spin_2_3_pt[i];
      }
    }
  }
}

void set_coarse_self_coupling_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, 
                                         vector_PRECISION *V, const int n, level_struct *l ) {
  
  int i, j, k, m, k1, k2, num_aggregates = l->is_PRECISION.num_agg,
    num_eig_vect = l->next_level->num_parent_eig_vect,
    aggregate_size = l->num_inner_lattice_sites*l->num_parent_eig_vect*2/num_aggregates,
    offset = l->num_parent_eig_vect,
    clover_site_size = (num_eig_vect*(2*num_eig_vect+1));
  buffer_PRECISION spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_PRECISION clover_pt, clover = l->next_level->op_PRECISION.clover;  
  
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for ( k=0; k<=n; k++ ) {
    k1 = (n*(n+1))/2+k; k2 = (n*(n+1))/2+k+(num_eig_vect*(num_eig_vect+1))/2;
  
    for ( j=0; j<num_aggregates; j++ ) {
      spin_0_1_pt = spin_0_1->vector_buffer + j*aggregate_size;
      spin_2_3_pt = spin_2_3->vector_buffer + j*aggregate_size;
      interpolation_data = V[k].vector_buffer + j*aggregate_size;
      clover_pt = clover + j*clover_site_size;
      
      for ( i=0; i<aggregate_size; ) {
        // A
        for ( m=0; m<offset; m++, i++ )
          clover_pt[ k1 ] += conj_PRECISION( interpolation_data[i] ) * spin_0_1_pt[i];
        // D
        for ( m=0; m<offset; m++, i++ )
          clover_pt[ k2 ] += conj_PRECISION( interpolation_data[i] ) * spin_2_3_pt[i];
      }
    }
  }
  
  for ( k=0; k<num_eig_vect; k++ ) {
    k1 = num_eig_vect*(num_eig_vect+1+n) + k;
  
    for ( j=0; j<num_aggregates; j++ ) {
      spin_0_1_pt = spin_0_1->vector_buffer + j*aggregate_size;
      spin_2_3_pt = spin_2_3->vector_buffer + j*aggregate_size;
      interpolation_data = V[k].vector_buffer + j*aggregate_size;
      clover_pt = clover + j*clover_site_size;
      
      for ( i=0; i<aggregate_size; ) {
        // B
        for ( m=0; m<offset; m++, i++ )
          clover_pt[ k1 ] += conj_PRECISION( interpolation_data[i] ) * spin_2_3_pt[i];
        i += offset;
      }
    }
  }
}


void set_coarse_neighbor_coupling_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, 
                                             vector_PRECISION *V, const int mu, const int n, level_struct *l ) {
  
  int i, i1, j, k, k1, k2, m, num_aggregates = l->is_PRECISION.num_agg,
      num_eig_vect = l->next_level->num_parent_eig_vect,
      offset = l->num_parent_eig_vect, nlsv = l->num_parent_eig_vect*2,
      D_link_size = num_eig_vect*num_eig_vect*4, *index_dir = l->is_PRECISION.agg_boundary_index[mu],
      aggregate_boundary_sites = l->is_PRECISION.agg_boundary_length[mu]/num_aggregates;
      
  buffer_PRECISION spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_PRECISION D_pt, D = l->next_level->op_PRECISION.D;
  
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  for ( k=0; k<num_eig_vect; k++ ) {
    k1 = n*num_eig_vect + k;
    k2 = (n+num_eig_vect)*num_eig_vect + k;
    i1 = 0;
    for ( j=0; j<num_aggregates; j++ ) {
      D_pt = D+(j*4+mu)*D_link_size;
      
      for ( i=0; i<aggregate_boundary_sites; i++ ) {
        spin_0_1_pt = spin_0_1->vector_buffer + nlsv*index_dir[i1];
        interpolation_data = V[k].vector_buffer + nlsv*index_dir[i1]; i1++;
        // A
        for ( m=0; m<offset; m++ )
          D_pt[ k1 ] += conj_PRECISION( interpolation_data[m] ) * spin_0_1_pt[m];
        // C
        for ( ; m<2*offset; m++ )
          D_pt[ k2 ] += conj_PRECISION( interpolation_data[m] ) * spin_0_1_pt[m];
      }
    }
    
    k1 = (n+2*num_eig_vect)*num_eig_vect + k;
    k2 = (n+3*num_eig_vect)*num_eig_vect + k;
    i1 = 0;
    for ( j=0; j<num_aggregates; j++ ) {
      D_pt = D+(j*4+mu)*D_link_size;
      
      for ( i=0; i<aggregate_boundary_sites; i++ ) {
        spin_2_3_pt = spin_2_3->vector_buffer + nlsv*index_dir[i1];
        interpolation_data = V[k].vector_buffer + nlsv*index_dir[i1]; i1++;
        // B
        for ( m=0; m<offset; m++ )
          D_pt[ k1 ] += conj_PRECISION( interpolation_data[m] ) * spin_2_3_pt[m];
        // D
        for ( ; m<2*offset; m++ )
          D_pt[ k2 ] += conj_PRECISION( interpolation_data[m] ) * spin_2_3_pt[m];
      }
    }
  }
}

void coarse_block_operator_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int n = s->num_block_sites, *length = s->dir_length, **index = s->index,
    *ind, *neighbor = s->op.neighbor_table, m = l->num_lattice_site_var, num_eig_vect = l->num_parent_eig_vect;
  vector_PRECISION lphi, leta;
  lphi.vector_buffer = phi->vector_buffer+start; leta.vector_buffer = eta->vector_buffer+start;
  vector_PRECISION leta1=leta, leta2=leta, lphi1=lphi, lphi2=lphi;

  // site-wise self coupling
  coarse_self_couplings_PRECISION( eta, phi, &(s->op), (start/m), (start/m)+n, l);

  // inner block couplings
  int hopp_size = 4 * SQUARE( num_eig_vect*2 );
  config_PRECISION D_pt, D = s->op.D + (start/m)*hopp_size;

  for ( int mu=0; mu<4; mu++ ) {
    ind = index[mu]; // mu direction
    for ( int i=0; i<length[mu]; i++ ) {
      int k = ind[i]; int j = neighbor[5*k+mu+1];
      D_pt = D + hopp_size*k + (hopp_size/4)*mu;
      leta1.vector_buffer = leta.vector_buffer+m*k;
      lphi1.vector_buffer = lphi.vector_buffer+m*j;
      coarse_hopp_PRECISION( &leta1, &lphi1, D_pt, l );

      leta2.vector_buffer = leta.vector_buffer+m*j;
      lphi2.vector_buffer = lphi.vector_buffer+m*k;
      coarse_daggered_hopp_PRECISION( &leta2, &lphi2, D_pt, l );
    }
  }
  END_UNTHREADED_FUNCTION(threading)
}


void coarse_aggregate_self_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
                                                schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, mu, index1, index2, length, *index_dir, *neighbor = s->op.neighbor_table,
      n = l->num_lattice_site_var, Dls = n*n, Dss = 4*n*n;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  config_PRECISION D_pt, D = s->op.D;
  
  vector_PRECISION_define( eta1, 0, 0, l->vector_size, l );
  vector_PRECISION_define( eta2, 0, 0, l->vector_size, l );  
  coarse_spinwise_self_couplings_PRECISION( eta1, eta2, phi, s->op.clover, l->inner_vector_size, l );
  
  for ( mu=0; mu<4; mu++ ) { // direction mu
    length = l->is_PRECISION.agg_length[mu]; index_dir = l->is_PRECISION.agg_index[mu];
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[5*index1+mu+1]; D_pt = D + Dss*index1 + Dls*mu;
      phi_pt.vector_buffer = phi->vector_buffer + n*index2; eta1_pt.vector_buffer = eta1->vector_buffer + n*index1; eta2_pt.vector_buffer = eta2->vector_buffer + n*index1;
      coarse_spinwise_n_hopp_PRECISION( &eta1_pt, &eta2_pt, &phi_pt, D_pt, l );
      phi_pt.vector_buffer = phi->vector_buffer + n*index1; eta1_pt.vector_buffer = eta1->vector_buffer + n*index2; eta2_pt.vector_buffer = eta2->vector_buffer + n*index2;
      coarse_spinwise_n_daggered_hopp_PRECISION( &eta1_pt, &eta2_pt, &phi_pt, D_pt, l );
    }
  }
}


void coarse_aggregate_neighbor_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
                                                    const int mu, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, index1, index2, length = l->is_PRECISION.agg_boundary_length[mu],
      *index_dir = l->is_PRECISION.agg_boundary_index[mu],
      *neighbor = l->is_PRECISION.agg_boundary_neighbor[mu],
      n = l->num_lattice_site_var, Dls = n*n, Dss = 4*n*n;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  config_PRECISION D_pt, D = s->op.D;
  
  vector_PRECISION_define( eta1, 0, 0, l->vector_size, l );
  vector_PRECISION_define( eta2, 0, 0, l->vector_size, l ); 
  
  // requires the positive boundaries of phi to be communicated befor
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i];
    index2 = neighbor[i];
    D_pt = D + Dss*index1 + Dls*mu;
    phi_pt.vector_buffer = phi->vector_buffer + n*index2; eta1_pt.vector_buffer = eta1->vector_buffer + n*index1; eta2_pt.vector_buffer = eta2->vector_buffer + n*index1;
    coarse_spinwise_hopp_PRECISION( &eta1_pt, &eta2_pt, &phi_pt, D_pt, l );
  }
}

void coarse_self_couplings_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
                                      operator_PRECISION_struct *op, int start, int end, level_struct *l ) {

  int num_eig_vect = l->num_parent_eig_vect, 
    vector_size = l->num_lattice_site_var,
    clover_size = (2*num_eig_vect*num_eig_vect+num_eig_vect), 
    block_size = (num_eig_vect*num_eig_vect+num_eig_vect);

  coarse_self_couplings_clover_PRECISION( eta+start*vector_size, phi+start*vector_size,
                                          op->clover+start*clover_size, (end-start)*vector_size, l );
#ifdef HAVE_TM // tm_term
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 )
    coarse_add_anti_block_diagonal_PRECISION( eta+start*vector_size, phi+start*vector_size, 
                                              op->tm_term+start*block_size, (end-start)*vector_size, l );
#endif
#ifdef HAVE_TM1p1 //eps_term
  if ( g.n_flavours == 2 &&
       ( op->epsbar != 0 || op->epsbar_ig5_odd_shift != 0 || op->epsbar_ig5_odd_shift != 0 ) )
    coarse_add_doublet_coupling_PRECISION( eta+start*vector_size, phi+start*vector_size, 
                                           op->epsbar_term+start*block_size, (end-start)*vector_size, l );
#endif

}

void coarse_aggregate_block_diagonal_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
                                                config_PRECISION block, level_struct *l ) {
  int length = l->inner_vector_size,
    num_eig_vect = l->num_parent_eig_vect,
    block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
  config_PRECISION block_pt = block;
  vector_PRECISION phi_pt=*phi, eta1_pt=*eta1, eta2_pt=*eta2, phi_end_pt;
  phi_end_pt.vector_buffer=phi->vector_buffer+length;
  // U(x) = [ A 0      , A=A*, D=D* 
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, columnwise
  // diagonal coupling
  while ( phi_pt.vector_buffer< phi_end_pt.vector_buffer ) {
    // A
    mvp_PRECISION( eta1_pt.vector_buffer, block_pt, phi_pt.vector_buffer, num_eig_vect );
    vector_PRECISION_define( &eta2_pt, _COMPLEX_PRECISION_ZERO, 0, num_eig_vect, l );
    block_pt += block_step_size; eta1_pt.vector_buffer += num_eig_vect; eta2_pt.vector_buffer += num_eig_vect; phi_pt.vector_buffer += num_eig_vect;
    // D
    vector_PRECISION_define( &eta1_pt, _COMPLEX_PRECISION_ZERO, 0, num_eig_vect, l );
    mvp_PRECISION( eta2_pt.vector_buffer, block_pt, phi_pt.vector_buffer, num_eig_vect );
    block_pt += block_step_size; eta1_pt.vector_buffer += num_eig_vect; eta2_pt.vector_buffer += num_eig_vect; phi_pt.vector_buffer += num_eig_vect;
  }
}



void coarse_spinwise_self_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
                                               config_PRECISION clover, int length, level_struct *l ) {
  
  int num_eig_vect = l->num_parent_eig_vect,
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,
      clover_step_size2 = SQUARE(num_eig_vect);
  config_PRECISION clover_pt = clover;
  buffer_PRECISION phi_pt=phi->vector_buffer, eta1_pt=eta1->vector_buffer, eta2_pt=eta2->vector_buffer+num_eig_vect, phi_end_pt=phi->vector_buffer+length;
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  while ( phi_pt < phi_end_pt ) {
    // A
    mvp_PRECISION( eta1_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size1; phi_pt += num_eig_vect; eta1_pt += num_eig_vect; 
    // D
    mvp_PRECISION( eta2_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size1; phi_pt -= num_eig_vect; eta2_pt -= num_eig_vect; 
    // C = -B*
    nmvh_PRECISION( eta1_pt, clover_pt, phi_pt, num_eig_vect );
    phi_pt += num_eig_vect; eta1_pt += num_eig_vect;
    // B
    mv_PRECISION( eta2_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size2; phi_pt += num_eig_vect; eta2_pt += 3*num_eig_vect;
  }
}

void coarse_operator_PRECISION_set_couplings( operator_PRECISION_struct *op, level_struct *l,
                                              struct Thread *threading ) {

  coarse_operator_PRECISION_set_neighbor_couplings( op, l, threading );
  coarse_operator_PRECISION_set_self_couplings( op, l, threading );

}

void coarse_operator_PRECISION_set_neighbor_couplings( operator_PRECISION_struct *op, level_struct *l,
                                                       struct Thread *threading ) {

}

void coarse_operator_PRECISION_set_self_couplings( operator_PRECISION_struct *op, level_struct *l, 
                                                   struct Thread *threading ) {
    
}

void coarse_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, int end, level_struct *l ) {
  
  int j, k=l->num_lattice_site_var/2;
  buffer_PRECISION eta_end, eta_pt, phi_pt;
  eta_end = eta->vector_buffer + end;
  phi_pt = phi->vector_buffer + start;
  eta_pt = eta->vector_buffer + start;
  
  if ( eta_pt != phi_pt ) {
    while ( eta_pt < eta_end ) {
      for ( j=0; j<k; j++ ) {
        *eta_pt = -(*phi_pt);
        eta_pt++; phi_pt++;
      }
      for ( j=0; j<k; j++ ) {
        *eta_pt = *phi_pt;
        eta_pt++; phi_pt++;
      }
    }
  } else {
    while ( eta_pt < eta_end ) {
      for ( j=0; j<k; j++ ) {
        *eta_pt = -(*eta_pt);
        eta_pt++;
      }
      eta_pt+=k;
    }
  }
}

void coarse_tau1_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, int end, level_struct *l ) {
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int j, k=l->num_lattice_site_var/4;
    buffer_PRECISION eta_end, phi_pt, eta_pt;
    
    eta_end = eta->vector_buffer + end;
    phi_pt = phi->vector_buffer + start;
    eta_pt = eta->vector_buffer + start;
    
    ASSERT( eta_pt != phi_pt );
    while ( eta_pt < eta_end ) {
      phi_pt += k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = -(*phi_pt);
        eta_pt++; phi_pt++;
      }
      phi_pt -= 2*k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = -(*phi_pt);
        eta_pt++; phi_pt++;
      }
      phi_pt += 2*k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = *phi_pt;
        eta_pt++; phi_pt++;
      }
      phi_pt -= 2*k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = *phi_pt;
        eta_pt++; phi_pt++;
      }
      phi_pt += k;
    }
  } else 
#endif
    {
      warning0("coarse_tau1_gamma5_PRECISION called with g.n_flavours != 2\n");
      coarse_gamma5_PRECISION( eta, phi, start, end, l );
    }
}

void apply_coarse_operator_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
                                      level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _SC, threading );
  int start;
  int end;
  compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &end, l, threading, 1);

  coarse_self_couplings_PRECISION( eta, phi, op, start, end, l);

  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );

  coarse_hopping_term_PRECISION( eta, phi, op, _FULL_SYSTEM, l, threading );

  PROF_PRECISION_STOP( _NC, 1, threading );
}

void g5D_apply_coarse_operator_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
                                          level_struct *l, struct Thread *threading ) {
  int start, end;
  compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
  apply_coarse_operator_PRECISION( eta, phi, op, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( eta, eta, start, end, l );
  SYNC_CORES(threading)
}


void apply_coarse_operator_dagger_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
                                             level_struct *l, struct Thread *threading ) {
  
  coarse_gamma5_PRECISION( &(l->vbuf_PRECISION[3]), phi, threading->start_index[l->depth], threading->end_index[l->depth], l );
  apply_coarse_operator_PRECISION( eta, &(l->vbuf_PRECISION[3]), op, l, threading );
  coarse_gamma5_PRECISION( eta, eta, threading->start_index[l->depth], threading->end_index[l->depth], l );
}


void coarse_operator_PRECISION_test_routine( level_struct *l, struct Thread *threading ) {

  if ( !l->idle ) {
    int ivs = l->inner_vector_size, civs = l->next_level->inner_vector_size;
    PRECISION diff = 0;
    vector_PRECISION vp[4], vc[3];
    
    for(int i=0; i<4; i++){
      vector_PRECISION_init( &vp[i] );
      vector_PRECISION_alloc( &vp[i], _ORDINARY, 1, l, threading );
    }

    for(int i=0; i<3; i++){                                                               
      vector_PRECISION_init( &vc[i] );                                                    
      vector_PRECISION_alloc( &vc[i], _ORDINARY, 1, l->next_level, threading );                       
    } 
    
    SYNC_MASTER_TO_ALL(threading)
    
    START_LOCKED_MASTER(threading)
#ifdef HAVE_TM1p1
    if(g.n_flavours == 1)
#endif
    {
      diff = global_inner_product_PRECISION( &(l->is_PRECISION.interpolation[0]), &(l->is_PRECISION.interpolation[1]), 0, ivs, l, no_threading )
        / global_norm_PRECISION( &(l->is_PRECISION.interpolation[0]), 0, ivs, l, no_threading );
      
      test0_PRECISION("depth: %d, correctness of block_gram_schmidt: %le\n", l->depth, cabs(diff) );
    }
    
    if ( !l->next_level->idle )
      vector_PRECISION_define_random( &vc[0], 0, civs, l->next_level );
    vector_PRECISION_distribute( &vc[1], &vc[0], l->next_level );
    vector_PRECISION_gather( &vc[2], &vc[1], l->next_level );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( &vc[1], &vc[0], &vc[2], 0, civs, l->next_level );
      diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
    }
    test0_PRECISION("depth: %d, correctness of gather( distribute( phi_c ) ) : %le\n", l->depth, diff );
        
    if ( !l->next_level->idle )
      vector_PRECISION_define_random( &vc[0], 0, civs, l->next_level );
    interpolate3_PRECISION( &vp[0], &vc[0], l, no_threading );
    restrict_PRECISION( &vc[1], &vp[0], l, no_threading );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( &vc[2], &vc[0], &vc[1], 0, civs, l->next_level );
      diff = global_norm_PRECISION( &vc[2], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
      test0_PRECISION("depth: %d, correctness of ( P* P - 1 ) phi_c: %le\n", l->depth, abs_PRECISION(diff) );
    }    
      
    END_LOCKED_MASTER(threading)
    if(threading->n_core>1) {
      interpolate3_PRECISION( &vp[0], &vc[0], l, threading );
      restrict_PRECISION( &vc[1], &vp[0], l, threading );
      START_LOCKED_MASTER(threading)
      if ( !l->next_level->idle ) {
        vector_PRECISION_minus( &vc[2], &vc[0], &vc[1], 0, civs, l->next_level );
        diff = global_norm_PRECISION( &vc[2], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
        test0_PRECISION("depth: %d, correctness of ( P* P - 1 ) phi_c with threading: %le\n", l->depth, diff );
      }
      END_LOCKED_MASTER(threading)
    }

    START_LOCKED_MASTER(threading)
    if (l->depth==0) 
      gamma5_PRECISION( &vp[1], &vp[0], l, no_threading );
    else
      coarse_gamma5_PRECISION( &vp[1], &vp[0], 0, ivs, l );
    restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
    coarse_gamma5_PRECISION( &vc[2], &vc[1], 0, civs, l->next_level );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( &vc[1], &vc[0], &vc[2], 0, civs, l->next_level );
      diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
      test0_PRECISION("depth: %d, correctness of ( g5_c P* g5 P - 1 ) phi_c: %le\n", l->depth, diff );
    }    
#ifdef HAVE_TM1p1
    if(g.n_flavours == 2) {
      if (l->depth==0) 
        tau1_gamma5_PRECISION( &vp[1], &vp[0], l, no_threading );
      else
        coarse_tau1_gamma5_PRECISION( &vp[1], &vp[0], 0, ivs, l );
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      coarse_tau1_gamma5_PRECISION( &vc[2], &vc[1], 0, civs, l->next_level );
      if ( !l->next_level->idle ) {
        vector_PRECISION_minus( &vc[1], &vc[0], &vc[2], 0, civs, l->next_level );
        diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
        test0_PRECISION("depth: %d, correctness of ( tau1 g5_c P* tau1 g5 P - 1 ) phi_c: %le\n", l->depth, diff );
      }    
    }
#endif
    END_LOCKED_MASTER(threading)

    START_LOCKED_MASTER(threading)
    vector_PRECISION_define( &vp[1], 0, 0, ivs, l );
    if (l->depth==0) 
      add_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.odd_proj, ivs );
    else
      coarse_add_block_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.odd_proj, ivs, l );
    restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
    
    vector_PRECISION_scale( &vc[1], &vc[1], -1.0, 0, civs, l->next_level );
    coarse_add_block_diagonal_PRECISION( &vc[1], &vc[0], l->next_level->s_PRECISION.op.odd_proj, civs, l->next_level );
    diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
    test0_PRECISION("depth: %d, correctness of ( P* 1odd P - 1odd_c ) phi_c: %le\n", l->depth, diff );
    END_LOCKED_MASTER(threading)  

#ifdef HAVE_TM
    START_LOCKED_MASTER(threading)
    if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) {
      vector_PRECISION_define( &vp[1], 0, 0, ivs, l );
      if (l->depth==0) 
        add_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.tm_term, ivs );
      else
        coarse_add_anti_block_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.tm_term, ivs, l );
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      
      vector_PRECISION_scale( &vc[1], &vc[1], -g.mu_factor[l->next_level->depth]/g.mu_factor[l->depth], 0, civs, l->next_level );
      coarse_add_anti_block_diagonal_PRECISION( &vc[1], &vc[0], l->next_level->s_PRECISION.op.tm_term, civs, l->next_level );
      diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
      test0_PRECISION("depth: %d, correctness of ( P* tm P - tm_c ) phi_c: %le\n", l->depth, diff );
    }
    END_LOCKED_MASTER(threading)  
#endif

#ifdef HAVE_TM1p1
    START_LOCKED_MASTER(threading)
    if ( g.n_flavours == 2 &&
         ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) ) {
      vector_PRECISION_define( &vp[1], 0, 0, ivs, l );
      if (l->depth==0) 
        apply_doublet_coupling_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.epsbar_term, ivs );
      else
        coarse_add_doublet_coupling_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.epsbar_term, ivs, l );
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      
      vector_PRECISION_scale( &vc[1], &vc[1], -g.epsbar_factor[l->next_level->depth]/g.epsbar_factor[l->depth], 0, civs, l->next_level );
      coarse_add_doublet_coupling_PRECISION( &vc[1], &vc[0], l->next_level->s_PRECISION.op.epsbar_term, civs, l->next_level );
      diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
      test0_PRECISION("depth: %d, correctness of ( P* eps P - eps_c ) phi_c: %le\n", l->depth, diff );
    }
    END_LOCKED_MASTER(threading)  
#endif

    if ( l->level > 0 ) {
      START_LOCKED_MASTER(threading)
      interpolate3_PRECISION( &vp[0], &vc[0], l, no_threading );

      apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );      
      
#ifdef HAVE_TM
      if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
        if (g.mu_factor[l->depth] != g.mu_factor[l->next_level->depth]) {  
          vector_PRECISION_scale( &vp[2], &vp[0], (g.mu_factor[l->next_level->depth]/g.mu_factor[l->depth])-1., 0, ivs, l );
          if(l->depth == 0)
            add_diagonal_PRECISION( &vp[1], &vp[2], l->p_PRECISION.op->tm_term, ivs );
          else
            coarse_add_anti_block_diagonal_PRECISION( &vp[1], &vp[2], l->p_PRECISION.op->tm_term, ivs, l );
        }
#endif      
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );

      if ( !l->next_level->idle ) {
        if ( l->level==1 && g.odd_even )
          coarse_odd_even_PRECISION_test( &vc[2], &vc[0], l->next_level, no_threading );
        else
          apply_operator_PRECISION( &vc[2], &vc[0], &(l->next_level->p_PRECISION), l->next_level, no_threading );
        
        vector_PRECISION_minus( &vc[2], &vc[1], &vc[2], 0, civs, l->next_level );
        diff = global_norm_PRECISION( &vc[2], 0, civs, l->next_level, no_threading ) /global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading );

        if ( l->level==1 && g.odd_even ) {
          test0_PRECISION("depth: %d, correctness of odd even preconditioned ( P* D P - D_c ) phi_c: %le\n", l->depth, diff );
        } else {
          test0_PRECISION("depth: %d, correctness of ( P* D P - D_c ) phi_c: %le\n", l->depth, diff );
        }
      }      
      END_LOCKED_MASTER(threading)

      if(threading->n_core>1) {
        if ( !l->next_level->idle ) {
          if ( l->level==1 && g.odd_even )
            coarse_odd_even_PRECISION_test( &vc[2], &vc[0], l->next_level, threading );
          else
            apply_operator_PRECISION( &vc[2], &vc[0], &(l->next_level->p_PRECISION), l->next_level, threading );
        }
        START_LOCKED_MASTER(threading)
        if ( !l->next_level->idle ) {
          vector_PRECISION_minus( &vc[2], &vc[1], &vc[2], 0, civs, l->next_level );
          diff = global_norm_PRECISION( &vc[2], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading );
          if ( l->level==1 && g.odd_even ) { //TODO: this test doesn't work without SSE!!
            test0_PRECISION("depth: %d, correctness of odd even preconditioned ( P* D P - D_c ) phi_c with D_c threaded: %le\n", l->depth, diff );
          } else {
            test0_PRECISION("depth: %d, correctness of ( P* D P - D_c ) phi_c with D_c threaded: %le\n", l->depth, diff );
          }
        }
        END_LOCKED_MASTER(threading)
        }
    }
    START_LOCKED_MASTER(threading)

    if ( l->level > 0 && l->depth > 0 && g.method == 3 && g.odd_even ) {
      vector_PRECISION_define_random( &vp[0], 0, ivs, l );
      block_to_oddeven_PRECISION( &vp[3], &vp[0], l, no_threading );
      coarse_diag_ee_PRECISION( &vp[2], &vp[3], &(l->oe_op_PRECISION), l, no_threading );
      coarse_diag_oo_PRECISION( &vp[2], &vp[3], &(l->oe_op_PRECISION), l, no_threading );
      coarse_hopping_term_PRECISION( &vp[2], &vp[3], &(l->oe_op_PRECISION), _FULL_SYSTEM, l, no_threading );
      oddeven_to_block_PRECISION( &vp[3], &vp[2], l, no_threading );
      apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );
      vector_PRECISION_minus( &vp[3], &vp[3], &vp[1], 0, ivs, l );
      diff = global_norm_PRECISION( &vp[3], 0, ivs, l, no_threading ) / global_norm_PRECISION( &vp[1], 0, ivs, l, no_threading );
      test0_PRECISION("depth: %d, correctness of odd even layout (smoother): %le\n", l->depth, diff );
     
      block_to_oddeven_PRECISION( &vp[3], &vp[0], l, no_threading );
      coarse_odd_even_PRECISION_test( &vp[2], &vp[3], l, no_threading );
      oddeven_to_block_PRECISION( &vp[3], &vp[2], l, no_threading );
      apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );
      vector_PRECISION_minus( &vp[3], &vp[3], &vp[1], 0, ivs, l );
      diff = global_norm_PRECISION( &vp[3], 0, ivs, l, no_threading ) / global_norm_PRECISION( &vp[1], 0, ivs, l, no_threading );
      test0_PRECISION("depth: %d, correctness of odd even preconditioned operator (smoother): %le\n", l->depth, diff );
    }
    
    for(int i=0; i<4; i++)
      vector_PRECISION_free( &vp[i], l, threading );

    for(int i=0; i<3; i++)
      vector_PRECISION_free( &vc[i], l->next_level, threading );
    
   END_LOCKED_MASTER(threading)
    
    if ( g.method != 6 && l->next_level->level > 0  && !l->next_level->idle ) {
      schwarz_PRECISION_mvm_testfun( &(l->next_level->s_PRECISION), l->next_level, threading );
    }

    if ( l->next_level->level > 0 && !l->next_level->idle )
      coarse_operator_PRECISION_test_routine( l->next_level, threading );
    
    SYNC_CORES(threading)
    SYNC_MASTER_TO_ALL(threading)
  }
}
