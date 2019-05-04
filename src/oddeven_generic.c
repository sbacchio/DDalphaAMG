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

void selfcoupling_cholesky_decomposition_PRECISION( const config_PRECISION output, config_double input ) {
  
/*********************************************************************************
* Performs a cholesky decomposition for a selfcoupling term.
* Input = [ A 0 ]   , A=A*, B=B*
*         [ 0 B ]
* Output ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major 
*                  (matlab notation). 
*********************************************************************************/  
  
  register int i, j, k;
  int n, offset[4] = {0,12,6,27};
  config_double in_pt;
  config_PRECISION out_pt = output;
  complex_PRECISION L[6][6], s;
  
  for ( n=0; n<2; n++ ) {
    // construct initial L = tril(A) for n=0, L = tril(B) for n=1, L row major
    in_pt = input+offset[2*n];
    for ( j=0; j<6; j++ ) {
      L[j][j] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    in_pt = input+offset[2*n+1];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        L[i][j] = (complex_PRECISION) conj_double(*in_pt); in_pt++;
      }
    }
    
    // calculate cholesky factor L
    for ( i=0; i<6; i++ ) {
      for ( j=0; j<=i; j++ ) {
        s = L[i][j];
        for ( k=0; k<j; k++ )
          s = s - L[i][k]*conj_PRECISION(L[j][k]);
        if ( i > j )
          L[i][j] = s / L[j][j];
        else if ( cabs_PRECISION(L[i][i]) > EPS_PRECISION ) 
          L[i][i] = csqrt_PRECISION(s);
        else {
          L[i][i] = csqrt_PRECISION(s);
        }
      }  
    }
    
    // output = tril(L) row major 
    for ( i=0; i<6; i++ ) {
      for ( j=0; j<=i; j++ ) {
        *out_pt = L[i][j]; out_pt++;
      }
    }
  }
}

#ifdef HAVE_TM
void selfcoupling_LU_decomposition_PRECISION( const config_PRECISION output, config_double input ) {

  /*********************************************************************************   
   * Performs a LU decomposition for a selfcoupling term.
   * Input = [ A 0 ]   , A=A*, B=B* (diagonals excluded)
   *         [ 0 B ]               
   * Input ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation).
   * Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without 
   *                  the diagonal of L which is equal to 1
   *********************************************************************************/

  register int i, j, k;
  int n;
  config_double in_pt;
  config_PRECISION out_pt;

  int offset[4] = {0,12,6,27};

  // construct initial L = A for n=0, L = B for n=1, L row major
  for ( n=0; n<2; n++ ) {

    out_pt = output + n*36;
    
    in_pt = input + offset[2*n];
    for ( j=0; j<6; j++ ) {
      out_pt[6*j+j] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    in_pt = input+offset[2*n+1];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        out_pt[6*j+i] = (complex_PRECISION) *in_pt;
        out_pt[6*i+j] = (complex_PRECISION) conj_double(*in_pt); in_pt++;
      }
    }
    
    // calculate LU
    for ( k=0; k<5; k++ ) {
      for ( i=k+1; i<6; i++ ) {
        out_pt[6*i+k] = out_pt[6*i+k]/out_pt[6*k+k]; // L: out(i,k) = out(i,k)/out(k,k)
        for ( j=k+1; j<6; j++ )
          out_pt[6*i+j] = out_pt[6*i+j]-out_pt[6*i+k]*out_pt[6*k+j]; // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
}
#endif

#ifdef HAVE_TM1p1
void selfcoupling_LU_doublet_decomposition_PRECISION( const config_PRECISION output, config_double input ) {

  /*********************************************************************************   
   * Performs a LU decomposition for a selfcoupling term.
   * Input = [ A 0 ]   , A=A*, B=B* (diagonals excluded)
   *         [ 0 B ]               
   * Input ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation).
   * Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without 
   *                  the diagonal of L which is equal to 1
   *********************************************************************************/

  register int i, j, k;
  int n;
  config_double in_pt;
  config_PRECISION out_pt = output;

  int offset[8] = {0,12,24,54,6,18,39,60};
  
  // construct initial L = A for n=0, L = B for n=1, L row major
  for ( n=0; n<2; n++ ) {

    out_pt = output + n*144;

    in_pt = input+offset[4*n];
    for ( j=0; j<6; j++ ) {
      out_pt[12*j+j] = (complex_PRECISION) *in_pt; in_pt++;
    }
    in_pt = input+offset[4*n+1];
    for ( j=0; j<6; j++ ) {
      out_pt[12*(j+6)+(j+6)] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    in_pt = input+offset[4*n+2];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        out_pt[12*(j+6)+(i+6)] = out_pt[12*j+i] = (complex_PRECISION) *in_pt;
        out_pt[12*(i+6)+(j+6)] = out_pt[12*i+j] = (complex_PRECISION) conj_double(*in_pt); in_pt++;
      }
    }
    
    in_pt = input+offset[4*n+3];
    for ( j=0; j<6; j++ ) {
      for ( i=0; i<6; i++ ) {
        out_pt[12*(j+6)+i] = out_pt[12*j+(i+6)] = _COMPLEX_PRECISION_ZERO;
      }
    }
    for ( j=0; j<6; j++ ) {
      out_pt[12*(j+6)+j] = out_pt[12*j+(j+6)] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    // calculate LU
    for ( k=0; k<11; k++ ) {
      for ( i=k+1; i<12; i++ ) {
        out_pt[12*i+k] = out_pt[12*i+k]/out_pt[12*k+k]; // L: out(i,k) = out(i,k)/out(k,k)
        for ( j=k+1; j<12; j++ )
          out_pt[12*i+j] = out_pt[12*i+j]-out_pt[12*i+k]*out_pt[12*k+j]; // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
}
#endif


static inline void LLH_perform_fwd_bwd_subs_PRECISION( vector_PRECISION *x, vector_PRECISION *b, config_PRECISION L,
						       int start, int end ) {

/*********************************************************************************
* Solves L*(L^H)*x = b for x, i.e., the clover coupling for a single lattice 
* site.
* - vector_PRECISION *b: Right hand side.
* - vector_PRECISION *x: Solution.
* - config_PRECISION L: Cholesky factor ( lower triangular matrix )
*********************************************************************************/
  
  register int id, i, j;
  int n;
  buffer_PRECISION x_pt = x->vector_buffer, b_pt = b->vector_buffer;
  x_pt += start; b_pt += start;

  L += (start/12)*42;
  
  for ( id=start; id<end; id+=12 ) {
    for ( n=0; n<2; n++ ) {
      // forward substitution with L
      for ( i=0; i<6; i++ ) {
        x_pt[i] = b_pt[i];
        for ( j=0; j<i; j++ ) {
          x_pt[i] = x_pt[i] - *L * x_pt[j]; L++;
        }
        x_pt[i] = x_pt[i] / *L; L++;
      }
      L -= 21;
      // backward substitution with L^H
      for ( i=5; i>=0; i-- ) {
        for ( j=i+1; j<6; j++ ) {
          x_pt[i] = x_pt[i] - conj_PRECISION(L[(j*(j+1))/2 + i]) * x_pt[j];
        }
        x_pt[i] = x_pt[i] / conj_PRECISION(L[(i*(i+1))/2 + i]);
      }
      x_pt+=6;
      b_pt+=6;
      L+=21;
    }
    x_pt+=12; b_pt+=12; L+=42;
  }
}


static inline void LU_perform_fwd_bwd_subs_PRECISION( vector_PRECISION *x, vector_PRECISION *b, config_PRECISION LU,
 						      int start, int end ) {

/*********************************************************************************
* Solves L*U*x = b for x, i.e., the clover coupling for a single lattice 
* site.
* - vector_PRECISION *b: Right hand side.
* - vector_PRECISION *x: Solution.
* - config_PRECISION L: Lower matrix from modified LU decomposition
* Note: U is given by u_{ii}=1, u_{ij}=l_{ji}* / l_{ii} 
*********************************************************************************/
  
  register int id, i, j, n;
  buffer_PRECISION x_pt = x->vector_buffer, b_pt = b->vector_buffer;
  x_pt += start; b_pt += start;

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    LU += (start/24)*288;
    for ( id=start; id<end; id+=24 ) {
      for ( n=0; n<2; n++ ) {
        // solve x = U^(-1) L^(-1) b
        // forward substitution with L
        for ( i=0; i<12; i++ ) {
          x_pt[i] = b_pt[i];
          for ( j=0; j<i; j++ ) {
            x_pt[i] = x_pt[i] - LU[i*12+j]*x_pt[j];
          }
        }
        // backward substitution with U
        for ( i=12-1; i>=0; i-- ) {
          for ( j=i+1; j<12; j++ ) {
            x_pt[i] = x_pt[i] - LU[i*12+j]*x_pt[j];
          }
          x_pt[i] = x_pt[i]/LU[i*(12+1)];
        }
        x_pt+=12;
        b_pt+=12;
        LU+=12*12;
      }
      x_pt+=24; b_pt+=24; LU+=288;
    }
  } else
#endif
    {
      LU += (start/12)*72;
      for ( id=start; id<end; id+=12 ) {
        for ( n=0; n<2; n++ ) {
          // solve x = U^(-1) L^(-1) b
          // forward substitution with L
          for ( i=0; i<6; i++ ) {
            x_pt[i] = b_pt[i];
            for ( j=0; j<i; j++ ) {
              x_pt[i] = x_pt[i] - LU[i*6+j]*x_pt[j];
            }
          }
          // backward substitution with U
          for ( i=6-1; i>=0; i-- ) {
            for ( j=i+1; j<6; j++ ) {
              x_pt[i] = x_pt[i] - LU[i*6+j]*x_pt[j];
            }
            x_pt[i] = x_pt[i]/LU[i*(6+1)];
          }
          x_pt+=6;
          b_pt+=6;
          LU+=6*6;
        }
      }
      x_pt+=12; b_pt+=12; LU+=72;
    }
}

static inline void LLH_multiply_PRECISION( vector_PRECISION *y, vector_PRECISION *x, config_PRECISION L,
					   int start, int end ) {

/*********************************************************************************
* Applies the clover coupling term to a vector, by multiplying L^H 
* and then L. 
* - vector_PRECISION *x: Input vector.
* - vector_PRECISION *y: Output vector.
* - config_PRECISION L: Cholesky factor ( lower triangular matrix )
*********************************************************************************/
  
  register int id, i, j;
  int n;
  complex_PRECISION z[6];
  buffer_PRECISION x_pt = x->vector_buffer, y_pt = y->vector_buffer;
  x_pt += start; y_pt += start;  
  
  L += (start/12)*42;
  
  for ( id=start; id<end; id+=12 ) {
    for ( n=0; n<2; n++ ) {
      // z = L^H x
      for ( j=0; j<6; j++ ) { // columns
        for ( i=0; i<j; i++ ) { // rows
          z[i] += conj_PRECISION(*L)*x_pt[j]; L++;
        }
        z[j] = conj_PRECISION(*L)*x_pt[j]; L++;
      }
      L-=21;
      // y = L*z;
      for ( i=0; i<6; i++ ) { // rows
        y_pt[i] = *L * z[0]; L++;
        for ( j=1; j<=i; j++ ) { // columns
          y_pt[i] += *L * z[j]; L++;
        }
      }
      x_pt+=6;
      y_pt+=6;
    }
    x_pt+=12; y_pt+=12; L+=42;
  }
}

static inline void LU_multiply_PRECISION( vector_PRECISION *y, vector_PRECISION *x, config_PRECISION LU,
					  int start, int end ) {

/*********************************************************************************
* Applies the clover coupling term to a vector, by multiplying L^H 
* and then L. 
* - vector_PRECISION *x: Input vector.
* - vector_PRECISION *y: Output vector.
* - config_PRECISION LU: LU decomposition
*********************************************************************************/

  register int id, i, j, n;
  buffer_PRECISION x_pt = x->vector_buffer, y_pt = y->vector_buffer;
  x_pt += start; y_pt += start;

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    LU += (start/24)*288;
    for ( id=start; id<end; id+=24 ) {
      for ( n=0; n<2; n++ ) {
	for ( i=0; i<12; i++ ) {
	  y_pt[i] = LU[i*(12+1)]*x_pt[i];
	  for ( j=i+1; j<12; j++ )
	    y_pt[i] += LU[i*12+j]*x_pt[j];
	}
	// multiplication with L
	for ( i=12-1; i>0; i-- )
	  for ( j=0; j<i; j++ )
	    y_pt[i] += LU[i*12+j]*y_pt[j];
	
	x_pt+=12;
	y_pt+=12;
	LU+=12*12;
      }
      x_pt+=24; y_pt+=24; LU+=288;
    }
  } else
#endif
    {
      LU += (start/12)*72;
      for ( id=start; id<end; id+=12 ) {
	for ( n=0; n<2; n++ ) {
	  for ( i=0; i<6; i++ ) {
	    y_pt[i] = LU[i*(6+1)]*x_pt[i];
	    for ( j=i+1; j<6; j++ )
	      y_pt[i] += LU[i*6+j]*x_pt[j];
	  }
	  // multiplication with L
	  for ( i=6-1; i>0; i-- )
	    for ( j=0; j<i; j++ )
	      y_pt[i] += LU[i*6+j]*y_pt[j];
	  
	  x_pt+=6;
	  y_pt+=6;
	  LU+=6*6;
	}
      }
      x_pt+=12; y_pt+=12; LU+=72;
    }
}


void diag_ee_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, 
                        level_struct *l, int start, int end ) {

/*********************************************************************************
* Applies the even-even block of the odd even decomposition to a vector.
* - vector_PRECISION *x: Input vector.
* - vector_PRECISION *y: Output vector.
*********************************************************************************/


#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    LU_multiply_PRECISION( y, x, op->clover_doublet_oo_inv, start, end);
  } else {
#endif
    if ( g.csw ) {
#if defined(HAVE_TM) 
      LU_multiply_PRECISION( y, x, op->clover, start, end);
#else
      LLH_multiply_PRECISION( y, x, op->clover, start, end );
#endif
    } else {
      config_PRECISION sc = op->clover + start;
      for ( int i=start; i<end;  ) {
        FOR12( y->vector_buffer[i] = x->vector_buffer[i]*sc[i]; i++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
}

// for debugging only
void diag_ee_inv_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op,
                            level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    int n1 = op->num_even_sites;
    LU_perform_fwd_bwd_subs_PRECISION( y, x, op->clover_doublet_oo_inv, 0, n1*24);
  } else {
#endif
    int i, n1 = op->num_even_sites;
    if ( g.csw ) {
      // diagonal blocks applied to the even sites
#ifndef HAVE_TM
      LLH_perform_fwd_bwd_subs_PRECISION( y, x, op->clover, 0, n1*12 );
#else
      LU_perform_fwd_bwd_subs_PRECISION( y, x, op->clover, 0, n1*12 );
#endif
    } else {
      config_PRECISION sc = op->clover;
      for ( i=0; i<n1*12;  ) {
        FOR12( y->vector_buffer[i] = x->vector_buffer[i]/sc[i]; i++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
  END_UNTHREADED_FUNCTION(threading)
}

// for debugging only
void diag_oo_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op,
    level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies the odd-odd block of the odd even decomposition to a vector.
* - vector_PRECISION *x: Input vector.
* - vector_PRECISION *y: Output vector.
*********************************************************************************/

  START_UNTHREADED_FUNCTION(threading)

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    int n1 = op->num_even_sites, n2 = op->num_odd_sites;
    LU_multiply_PRECISION( y, x, op->clover_doublet_oo_inv, n1*24, (n1+n2)*24 );
  } else {
#endif
    int i, n1 = op->num_even_sites, n2 = op->num_odd_sites;
    // diagonal blocks applied to the odd sites
    if ( g.csw ) {
#ifndef HAVE_TM
      LLH_multiply_PRECISION( y, x, op->clover, n1*12, (n1+n2)*12 );
#else
      LU_multiply_PRECISION( y, x, op->clover, n1*12, (n1+n2)*12 );
#endif
    } else {
      config_PRECISION sc = op->clover + n1*12;
      for ( i=n1*12; i<(n1+n2)*12;  ) {
        FOR12( y->vector_buffer[i] = x->vector_buffer[i]*sc[i]; i++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif

  END_UNTHREADED_FUNCTION(threading)
}


void diag_oo_inv_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op,
                            level_struct *l, int start, int end ) {

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    // inverse diagonal blocks applied to the odd sites
    LU_perform_fwd_bwd_subs_PRECISION( y, x, op->clover_doublet_oo_inv, start, end );
  } else {
#endif
    // inverse diagonal blocks applied to the odd sites
    if ( g.csw ) {
#if defined(HAVE_TM)
      LU_perform_fwd_bwd_subs_PRECISION( y, x, op->clover, start, end );
#else
      LLH_perform_fwd_bwd_subs_PRECISION( y, x, op->clover, start, end );
#endif
    } else {
      config_PRECISION sc = op->clover + start;
      for ( int i=start; i<end; ) {
        FOR12( y->vector_buffer[i] = x->vector_buffer[i]/sc[i]; i++;)
      }
    }
#ifdef HAVE_TM1p1
  }
#endif
}


void oddeven_setup_PRECISION( operator_double_struct *in, level_struct *l ) {

/*********************************************************************************
* Reorder data layouts and index tables to allow for odd even preconditioning.
*********************************************************************************/ 
  
  int j, k, k_e, k_o, n=l->num_inner_lattice_sites, oe_offset=0, mu, nu,
    sc_size = g.csw ? 42:12, lu_dec_size = 42, bs, **bt = NULL,
      *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  config_double sc_in = in->clover, nc_in = in->D;
  config_PRECISION Aee = NULL, Aoo = NULL;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  op->m0 = in->m0;

#ifdef HAVE_TM
  op->mu = in->mu;
  op->mu_even_shift = in->mu_even_shift;
  op->mu_odd_shift = in->mu_odd_shift;

  lu_dec_size = 72;
  config_double tm_term_in = in->tm_term;
#endif
  
  for ( mu=0; mu<4; mu++ ) {
    le[mu] = l->local_lattice[mu];
    N[mu] = le[mu]+1;
    op->table_dim[mu] = N[mu];
  }
  
  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;
  
  // estimate site numbers
  op->num_even_sites = 0;
  op->num_odd_sites = 0;
  op->oe_offset = oe_offset;
  
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            op->num_odd_sites++;
          } else {
            op->num_even_sites++;
          }
        }
  
  // re-order clover term (i.e., self coupling)
  if ( g.csw ) {
    MALLOC( op->clover, complex_PRECISION, lu_dec_size*n );
    Aee = op->clover;
    Aoo = op->clover + op->num_even_sites*lu_dec_size;
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              // odd sites
#ifndef HAVE_TM
              selfcoupling_cholesky_decomposition_PRECISION( Aoo, sc_in );
#else
              complex_double buffer[42];
              for(int i=0; i<42; i++)
                buffer[i] = sc_in[i];
              for(int i=0; i<12; i++)
                buffer[i] += tm_term_in[i];
              selfcoupling_LU_decomposition_PRECISION( Aoo, buffer );
#endif
              Aoo += lu_dec_size;
            } else {
              // even sites
#ifndef HAVE_TM
              selfcoupling_cholesky_decomposition_PRECISION( Aee, sc_in );
#else
              complex_double buffer[42];
              for(int i=0; i<42; i++)
                buffer[i]=(complex_double) sc_in[i];
              for(int i=0; i<12; i++)
                buffer[i]+=(complex_double) tm_term_in[i];
              selfcoupling_LU_decomposition_PRECISION( Aee, buffer );
#endif
              Aee += lu_dec_size;
            }
            sc_in += sc_size;
#ifdef HAVE_TM
            tm_term_in += 12;
#endif
          }
  } else {
    MALLOC( op->clover, complex_PRECISION, 12*n );
    Aee = op->clover;
    Aoo = op->clover + op->num_even_sites*12;
    
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              // odd sites
              FOR12( *Aoo = *sc_in; Aoo++; sc_in++; )
            } else {
              // even sites
              FOR12( *Aee = *sc_in; Aee++; sc_in++; )
            }
          }
  }

#ifdef HAVE_TM1p1

  int lu_doublet_dec_size = 288;
  config_double eps_term_in = in->epsbar_term;
  sc_in = in->clover; 
#ifdef HAVE_TM
  tm_term_in = in->tm_term;
#endif
  op->epsbar = in->epsbar;
  op->epsbar_ig5_even_shift = in->epsbar_ig5_even_shift;
  op->epsbar_ig5_odd_shift = in->epsbar_ig5_odd_shift;
    
  // re-order clover term (i.e., self coupling)
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, lu_doublet_dec_size*n );
  Aee = op->clover_doublet_oo_inv;
  Aoo = op->clover_doublet_oo_inv + op->num_even_sites*lu_doublet_dec_size;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            // odd sites
            complex_double buffer[66];
            if ( g.csw ) {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = (complex_double) sc_in[i];
            } else {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = _COMPLEX_double_ZERO;
            }              
            for(int i=0; i<12; i++) //54-65
              buffer[i+54] = (complex_double) eps_term_in[i];
#ifdef HAVE_TM
            if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
              for(int i=0; i<12; i++) {
                buffer[i] += (complex_double) tm_term_in[i];
                buffer[i+12] -= (complex_double) tm_term_in[i];
              }
#endif
            selfcoupling_LU_doublet_decomposition_PRECISION( Aoo, buffer );
            Aoo += lu_doublet_dec_size;
          } else {
            // even sites
            complex_double buffer[66];
            if ( g.csw ) {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = (complex_double) sc_in[i];
            } else {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = _COMPLEX_double_ZERO;
            }              
            for(int i=0; i<12; i++) //54-65
              buffer[i+54] = (complex_double) eps_term_in[i];
#ifdef HAVE_TM
            if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
              for(int i=0; i<12; i++) {
                buffer[i] += (complex_double) tm_term_in[i];
                buffer[i+12] -= (complex_double) tm_term_in[i];
              }
#endif
            selfcoupling_LU_doublet_decomposition_PRECISION( Aee, buffer );
            Aee += lu_doublet_dec_size;
          }
          sc_in += sc_size;
          eps_term_in += 12;
#ifdef HAVE_TM
          tm_term_in += 12;
#endif
        }
#endif
  
  // re-order hopping term (i.e., nearest neighbor coupling)
  MALLOC( op->D, complex_PRECISION, 36*n );
  
  k=0; k_e=0; k_o=0;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            for ( j=0; j<36; j++)
              (op->D+(k_e+op->num_even_sites)*36)[j] = (complex_PRECISION) (nc_in+k*36)[j];
            k_e++;
          } else {
            for ( j=0; j<36; j++)
              (op->D+k_o*36)[j] = (complex_PRECISION) (nc_in+k*36)[j];
            k_o++;
          }
          k++;
        }
       
  // define data layout
  MALLOC( op->index_table, int, N[T]*N[Z]*N[Y]*N[X] );
  eot = op->index_table;
  
  define_eot( eot, N, l );
  
  // neighbor table, translation table
  MALLOC( op->neighbor_table, int, 5*N[T]*N[Z]*N[Y]*N[X] );
  MALLOC( op->backward_neighbor_table, int, 5*N[T]*N[Z]*N[Y]*N[X] );
  MALLOC( op->translation_table, int, le[T]*le[Z]*le[Y]*le[X] );
  nt = op->neighbor_table;
  tt = op->translation_table;
  
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );
  
  // boundary table
  for ( mu=0; mu<4; mu++ ) {
    bs = 1;
    le[mu] = 1;
    for ( nu=0; nu<4; nu++ )
      bs *= le[nu];
    
    MALLOC( op->c.boundary_table[2*mu], int, bs );
    op->c.boundary_table[2*mu+1] = op->c.boundary_table[2*mu];
    
    le[mu] = l->local_lattice[mu];
  }
  
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );
  
  j = (l->num_lattice_site_var/2)*l->num_lattice_sites;
#ifdef HAVE_TM1p1
  j *= 2;
#endif
  MALLOC( op->prnT, complex_PRECISION, j*8 );
  op->prnZ = op->prnT + j; op->prnY = op->prnZ + j; op->prnX = op->prnY + j;
  op->prpT = op->prnX + j; op->prpZ = op->prpT + j; op->prpY = op->prpZ + j; op->prpX = op->prpY + j;  
  MALLOC( op->buffer, vector_PRECISION, 2 );
  for(int i=0; i<2; i++ ){
    vector_PRECISION_init( &(op->buffer[i]) );
#ifdef HAVE_TM1p1
    vector_PRECISION_alloc( &(op->buffer[i]), _ORDINARY, 2, l, no_threading );
#else
    vector_PRECISION_alloc( &(op->buffer[i]), _ORDINARY, 1, l, no_threading );
#endif
  }
  ghost_alloc_PRECISION( 0, &(op->c), l );
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;
  l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
}


void oddeven_free_PRECISION( level_struct *l ) {
  
  int mu, nu, nc_size = 9, lu_dec_size = 42,
      *ll = l->local_lattice, n = l->num_inner_lattice_sites, bs;
#ifdef HAVE_TM
  lu_dec_size = 72;
#endif
      
  ghost_free_PRECISION( &(l->oe_op_PRECISION.c), l );
  FREE( l->oe_op_PRECISION.D, complex_PRECISION, 4*nc_size*n );
  if ( g.csw )
    FREE( l->oe_op_PRECISION.clover, complex_PRECISION, lu_dec_size*n );
  else
    FREE( l->oe_op_PRECISION.clover, complex_PRECISION, 12*n );
  FREE( l->oe_op_PRECISION.index_table, int, (ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.neighbor_table, int, 5*(ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.backward_neighbor_table, int, 5*(ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.translation_table, int, ll[T]*ll[Z]*ll[Y]*ll[X] );
  
  for ( mu=0; mu<4; mu++ ) {
    bs = 1;
    for ( nu=0; nu<4; nu++ )
      if ( mu != nu )
        bs *= ll[nu];
    
    FREE( l->oe_op_PRECISION.c.boundary_table[2*mu], int, bs );
    l->oe_op_PRECISION.c.boundary_table[2*mu+1] = NULL;
  }

  for(int i=0; i<2; i++ ){
#ifdef HAVE_TM1p1
    vector_PRECISION_free( &(l->oe_op_PRECISION.buffer[i]), l, no_threading );
#else
    vector_PRECISION_free( &(l->oe_op_PRECISION.buffer[i]), l, no_threading );
#endif
  }
  FREE( l->oe_op_PRECISION.buffer, vector_PRECISION, 2 );
#ifdef HAVE_TM1p1
  FREE( l->oe_op_PRECISION.prnT, complex_PRECISION, 2*(l->num_lattice_site_var/2)*l->num_lattice_sites*8 );
  FREE( l->oe_op_PRECISION.clover_doublet_oo_inv, complex_PRECISION, 288*n );
#else
  FREE( l->oe_op_PRECISION.prnT, complex_PRECISION, (l->num_lattice_site_var/2)*l->num_lattice_sites*8 );
#endif
}


void oddeven_to_serial_PRECISION( vector_double *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Translates a vector from an odd even PRECISION precision layout to a serial 
* double precision layout.
*********************************************************************************/ 
  
  int i, j, k,
      nsv = l->num_lattice_site_var, *tt = l->oe_op_PRECISION.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt[i];
    for ( j=0; j<nsv; j++ ) {
      out->vector_buffer[i*nsv+j] = (complex_double) in->vector_buffer[k*nsv+j];
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}


void serial_to_oddeven_PRECISION( vector_PRECISION *out, vector_double *in, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Translates a vector from a serial double precision layout to an odd even
* PRECISION precision layout.
*********************************************************************************/ 
  
  int i, j, k,
      nsv = l->num_lattice_site_var, *tt = l->oe_op_PRECISION.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt[i];
    for ( j=0; j<nsv; j++ ) {
      out->vector_buffer[k*nsv+j] = (complex_PRECISION) in->vector_buffer[i*nsv+j];
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}


void oddeven_to_block_PRECISION( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {

  int i, j, k, m,
      nsv = l->num_lattice_site_var, *tt_oe = l->oe_op_PRECISION.translation_table,
      *tt_b = l->s_PRECISION.op.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt_oe[i]; m = tt_b[i];
    for ( j=0; j<nsv; j++ ) {
      out->vector_buffer[m*nsv+j] = in->vector_buffer[k*nsv+j];
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}


void block_to_oddeven_PRECISION( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {

  int i, j, k, m,
      nsv = l->num_lattice_site_var, *tt_oe = l->oe_op_PRECISION.translation_table,
      *tt_b = l->s_PRECISION.op.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt_oe[i]; m = tt_b[i];
    for ( j=0; j<nsv; j++ ) {
      out->vector_buffer[k*nsv+j] = in->vector_buffer[m*nsv+j];
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}

void hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
                             const int amount, level_struct *l, struct Thread *threading ) {

  int start_even, end_even, start_odd, end_odd, n = l->num_inner_lattice_sites,
    *neighbor = op->neighbor_table, start=0, plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  
  SYNC_CORES(threading)  
    
  if ( amount == _EVEN_SITES || amount == _ODD_SITES ) {
    compute_core_start_end_custom(0, op->num_even_sites, &start_even, &end_even, l, threading, 1 );
    compute_core_start_end_custom(op->num_even_sites, op->num_even_sites+op->num_odd_sites, &start_odd, &end_odd, l, threading, 1 );
  } else {
    compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &n, l, threading, 1 );
  }
  
  SYNC_CORES(threading)  
  
  if ( amount == _EVEN_SITES ) {
    start = start_odd, n = end_odd;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    start = start_even, n = end_even;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  int i, *nb_pt;
  buffer_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;
  
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    // project in negative directions
    complex_PRECISION pbuf[12];
    for ( i=12*start, phi_pt=phi->vector_buffer+24*start; i<12*n; i+=12, phi_pt+=24 ) {
      dprp_T_PRECISION( op->prnT+i, phi_pt );
      dprp_Z_PRECISION( op->prnZ+i, phi_pt );
      dprp_Y_PRECISION( op->prnY+i, phi_pt );
      dprp_X_PRECISION( op->prnX+i, phi_pt );
    }
    // start communication in negative direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
    END_LOCKED_MASTER(threading) 
    // project plus dir and multiply with U dagger
    for ( phi_pt=phi->vector_buffer+24*start, end_pt=phi->vector_buffer+24*n, D_pt = op->D+36*start, nb_pt=neighbor+4*start; phi_pt<end_pt; phi_pt+=24 ) {
      // T dir
      i = 12*(*nb_pt); nb_pt++;
      dprn_T_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpT+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpT+i+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpT+i+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpT+i+9, D_pt, pbuf+9 ); D_pt += 9;
      // Z dir
      i = 12*(*nb_pt); nb_pt++;
      dprn_Z_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpZ+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpZ+i+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpZ+i+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpZ+i+9, D_pt, pbuf+9 ); D_pt += 9;
      // Y dir
      i = 12*(*nb_pt); nb_pt++;
      dprn_Y_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpY+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpY+i+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpY+i+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpY+i+9, D_pt, pbuf+9 ); D_pt += 9;
      // X dir
      i = 12*(*nb_pt); nb_pt++;
      dprn_X_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpX+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpX+i+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpX+i+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpX+i+9, D_pt, pbuf+9 ); D_pt += 9;
    }
    if ( amount == _EVEN_SITES ) {
      start = start_even, n = end_even;
    } else if ( amount == _ODD_SITES ) {
      start = start_odd, n = end_odd;
    }  
    // start communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
    // wait for communication in negative direction
    ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
    ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
    ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
    ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
    END_LOCKED_MASTER(threading) 
    // multiply with U and lift up minus dir
    for ( eta_pt=eta->vector_buffer+24*start, end_pt=eta->vector_buffer+24*n, D_pt = op->D+36*start, nb_pt=neighbor+4*start; eta_pt<end_pt; eta_pt+=24 ) {
      // T dir
      i = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnT+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnT+i+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnT+i+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnT+i+9 );
      dpbp_su3_T_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Z dir
      i = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnZ+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnZ+i+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnZ+i+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnZ+i+9 );
      dpbp_su3_Z_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Y dir
      i = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnY+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnY+i+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnY+i+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnY+i+9 );
      dpbp_su3_Y_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // X dir
      i = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnX+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnX+i+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnX+i+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnX+i+9 );
      dpbp_su3_X_PRECISION( pbuf, eta_pt ); D_pt += 9;
    }
    // wait for communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
    ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
    ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
    ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
    END_LOCKED_MASTER(threading) 
    // lift up plus dir
    for ( i=12*start, eta_pt=eta->vector_buffer+24*start; i<12*n; i+=12, eta_pt+=24 ) {
      dpbn_su3_T_PRECISION( op->prpT+i, eta_pt );
      dpbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
      dpbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
      dpbn_su3_X_PRECISION( op->prpX+i, eta_pt );
    }
  } else {
#endif
    // project in negative directions
    complex_PRECISION pbuf[6];
    for ( i=6*start, phi_pt=phi->vector_buffer+12*start; i<6*n; i+=6, phi_pt+=12 ) {
      prp_T_PRECISION( op->prnT+i, phi_pt );
      prp_Z_PRECISION( op->prnZ+i, phi_pt );
      prp_Y_PRECISION( op->prnY+i, phi_pt );
      prp_X_PRECISION( op->prnX+i, phi_pt );
    }
    // start communication in negative direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
    END_LOCKED_MASTER(threading) 
    // project plus dir and multiply with U dagger
    for ( phi_pt=phi->vector_buffer+12*start, end_pt=phi->vector_buffer+12*n, D_pt = op->D+36*start, nb_pt=neighbor+4*start; phi_pt<end_pt; phi_pt+=12 ) {
      // T dir
      i = 6*(*nb_pt); nb_pt++;
      prn_T_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpT+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpT+i+3, D_pt, pbuf+3 ); D_pt += 9;
      // Z dir
      i = 6*(*nb_pt); nb_pt++;
      prn_Z_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpZ+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpZ+i+3, D_pt, pbuf+3 ); D_pt += 9;
      // Y dir
      i = 6*(*nb_pt); nb_pt++;
      prn_Y_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpY+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpY+i+3, D_pt, pbuf+3 ); D_pt += 9;
      // X dir
      i = 6*(*nb_pt); nb_pt++;
      prn_X_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpX+i, D_pt, pbuf );
      mvmh_PRECISION( op->prpX+i+3, D_pt, pbuf+3 ); D_pt += 9;
    }
    if ( amount == _EVEN_SITES ) {
      start = start_even, n = end_even;
    } else if ( amount == _ODD_SITES ) {
      start = start_odd, n = end_odd;
    }  
    // start communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
    ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
    // wait for communication in negative direction
    ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
    ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
    ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
    ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
    END_LOCKED_MASTER(threading) 
    // multiply with U and lift up minus dir
    for ( eta_pt=eta->vector_buffer+12*start, end_pt=eta->vector_buffer+12*n, D_pt = op->D+36*start, nb_pt=neighbor+4*start; eta_pt<end_pt; eta_pt+=12 ) {
      // T dir
      i = 6*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnT+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnT+i+3 );
      pbp_su3_T_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Z dir
      i = 6*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnZ+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnZ+i+3 );
      pbp_su3_Z_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Y dir
      i = 6*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnY+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnY+i+3 );
      pbp_su3_Y_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // X dir
      i = 6*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnX+i );
      mvm_PRECISION( pbuf+3, D_pt, op->prnX+i+3 );
      pbp_su3_X_PRECISION( pbuf, eta_pt ); D_pt += 9;
    }
    // wait for communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
    ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
    ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
    ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
    END_LOCKED_MASTER(threading) 
    // lift up plus dir
    for ( i=6*start, eta_pt=eta->vector_buffer+12*start; i<6*n; i+=6, eta_pt+=12 ) {
      pbn_su3_T_PRECISION( op->prpT+i, eta_pt );
      pbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
      pbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
      pbn_su3_X_PRECISION( op->prpX+i, eta_pt );
    }
#ifdef HAVE_TM1p1
  }
#endif

  SYNC_CORES(threading)
}

void apply_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
    level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies the Schur complement to a vector.
*********************************************************************************/

  // start and end indices for vector functions depending on thread
  int start_even, end_even, start_odd, end_odd;
  
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );
  
  vector_PRECISION *tmp = op->buffer;
  
  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp[0], 0, start_odd, end_odd, l );
  vector_PRECISION_define( &tmp[0], 0, start_even, end_even, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  
  PROF_PRECISION_START( _SC, threading );
  diag_ee_PRECISION( out, in, op, l, start_even, end_even );
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SC, 1, threading );
  hopping_term_PRECISION( &tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &tmp[1], &tmp[0], op, l, start_odd, end_odd );
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &tmp[0], &tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  vector_PRECISION_minus( out, out, &tmp[0], start_even, end_even, l );
}


void solve_oddeven_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  compute_core_start_end(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading);

  vector_PRECISION tmp = op->buffer[0];
  
  // odd to even
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &tmp, &(p->b), op, l, start, end );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
  vector_PRECISION_scale( &tmp, &tmp, -1, start, end, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &(p->b), &tmp, op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  
  if ( g.method == 4 )
    fgmres_PRECISION( p, l, threading );
  else if ( g.method == 5 )
    bicgstab_PRECISION( p, l, threading );
  diag_oo_inv_PRECISION( &(p->x), &(p->b), op, l, start, end );
  
  // even to odd
  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp, 0, start, end, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &tmp, &(p->x), op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &(p->b), &tmp, op, l, start, end );
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
  vector_PRECISION_minus( &(p->x), &(p->x), &(p->b), start, end, l );
  SYNC_CORES(threading)
}


void g5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, int end, level_struct *l ) {
  if ( eta->vector_buffer != phi->vector_buffer ) { 
    for( int i = start; i < end; ) {
      FOR6( eta->vector_buffer[i] = -phi->vector_buffer[i]; i++; )
      FOR6( eta->vector_buffer[i] = phi->vector_buffer[i]; i++; )
    }
  } else {
    for ( int i = start; i < end; ) {
      FOR6( eta->vector_buffer[i] = phi->vector_buffer[i]; i++; )
      i+=6;
    }
  }
}


void minus_g5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, int end, level_struct *l ) {
  if ( eta->vector_buffer != phi->vector_buffer ) { 
    for ( int i = start; i < end; ) {
      FOR6( eta->vector_buffer[i] = phi->vector_buffer[i]; i++; )
      FOR6( eta->vector_buffer[i] = -phi->vector_buffer[i]; i++; )
    }
  } else {
    for ( int i = start; i < end;  ) {
      i+=6;
      FOR6( eta->vector_buffer[i] = -phi->vector_buffer[i]; i++; )
    }
  }
}


void g5D_apply_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies the Schur complement to a vector.
*********************************************************************************/

  // start and end indices for vector functions depending on thread
  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );
  
  vector_PRECISION *tmp = op->buffer;

 // vector_PRECISION **tmp;
 // *tmp->vector_buffer = op->buffer->vector_buffer;
  

  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp[0], 0, start_odd, end_odd, l );
  vector_PRECISION_define( &tmp[0], 0, start_even, end_even, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  
  PROF_PRECISION_START( _SC, threading );
  diag_ee_PRECISION( out, in, op, l, start_even, end_even );
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SC, 1, threading );
  hopping_term_PRECISION( &tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &tmp[1], &tmp[0], op, l, start_odd, end_odd );
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &tmp[0], &tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  vector_PRECISION_minus( out, out, &tmp[0], start_even, end_even, l );
  SYNC_CORES(threading)
  g5_PRECISION( out, out, start_even, end_even, l );
//   g5_PRECISION( out, out, start_odd, end_odd, l );
  SYNC_CORES(threading)
}


void g5D_solve_oddeven_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );

  vector_PRECISION tmp = op->buffer[0];
  
  // odd to even
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &tmp, &(p->b), op, l, start_odd, end_odd );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
//   g5_PRECISION( tmp, tmp, start_odd, end_odd, l );
//   vector_PRECISION_scale( tmp, tmp, -1, start_odd, end_odd, l );
  minus_g5_PRECISION( &tmp, &tmp, start_odd, end_odd, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  vector_PRECISION_define( &(p->x), 0, start_even, end_even, l );
  hopping_term_PRECISION( &(p->x), &tmp, op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  SYNC_CORES(threading)
  g5_PRECISION( &(p->x), &(p->x), start_even, end_even, l );
  vector_PRECISION_plus( &(p->b), &(p->b), &(p->x), start_even, end_even, l );
  SYNC_CORES(threading)
  
  ASSERT( g.method == 6 );
  fgmres_PRECISION( p, l, threading );
  diag_oo_inv_PRECISION( &(p->x), &(p->b), op, l, start_odd, end_odd );
  g5_PRECISION( &(p->x), &(p->x), start_odd, end_odd, l );
  
  // even to odd
  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp, 0, start_odd, end_odd, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &tmp, &(p->x), op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &(p->b), &tmp, op, l, start_odd, end_odd );
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
  vector_PRECISION_minus( &(p->x), &(p->x), &(p->b), start_odd, end_odd, l );
  SYNC_CORES(threading)
}

// ----- block odd even -----------------------------------------------------------

void schwarz_PRECISION_oddeven_setup( schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int mu, i, d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, agg_split[4], block_split[4], block_size[4];
  operator_PRECISION_struct *op = &(s->op);
  int n1 = s->num_block_even_sites;
#ifdef HAVE_TM
  config_PRECISION tm_term_pt = op->tm_term;
#endif

  for ( mu=0; mu<4; mu++ ) {
    agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];
    block_split[mu] = l->coarsening[mu]/l->block_lattice[mu];
    block_size[mu] = l->block_lattice[mu];
  }
  
  if ( g.csw ) {
    config_PRECISION clover_pt = op->clover, clover_oo_inv_pt = op->clover_oo_inv;
    complex_double buffer[42];
    int cs = 42;
    for ( d0=0; d0<agg_split[T]; d0++ )
      for ( c0=0; c0<agg_split[Z]; c0++ )
        for ( b0=0; b0<agg_split[Y]; b0++ )
          for ( a0=0; a0<agg_split[X]; a0++ )
            
            for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
              for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
                for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                  for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {

                    // skipping even sites
                    clover_pt += n1*cs;
#ifdef HAVE_TM
                    tm_term_pt += n1*12;
#endif
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                              for ( i=0; i<42; i++ )
                                buffer[i] = (complex_double)clover_pt[i];
#ifdef HAVE_TM
                              if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
                                for ( i=0; i<12; i++ )
                                  buffer[i] += (complex_double)tm_term_pt[i];
                              tm_term_pt += 12;
                              selfcoupling_LU_decomposition_PRECISION( clover_oo_inv_pt, buffer );
                              clover_oo_inv_pt += 72;
#else  
                              selfcoupling_cholesky_decomposition_PRECISION( clover_oo_inv_pt, buffer );
                              clover_oo_inv_pt += 42;
#endif
                              clover_pt += cs;

                            }
                          }
                  }
  }

#ifdef HAVE_TM1p1
  complex_double buffer[66];
  config_PRECISION clover_oo_inv_pt = op->clover_doublet_oo_inv, clover_pt = op->clover;
  int cs = g.csw ? 42:12;
  config_PRECISION eps_term_pt = op->epsbar_term;
#ifdef HAVE_TM
  tm_term_pt = op->tm_term;
#endif
  
  for ( d0=0; d0<agg_split[T]; d0++ )
    for ( c0=0; c0<agg_split[Z]; c0++ )
      for ( b0=0; b0<agg_split[Y]; b0++ )
        for ( a0=0; a0<agg_split[X]; a0++ )
          
          for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
            for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
              for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {

                  // skipping even sites
                  clover_pt += n1*cs;
                  eps_term_pt += n1*12;
#ifdef HAVE_TM
                  tm_term_pt += n1*12;
#endif
                  for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                    for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                      for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                        for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                          if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                               (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                            if ( g.csw ) {
                              for( i=0; i<12; i++ ) //0-23
                                buffer[i+12] = buffer[i] = (complex_double) clover_pt[i];
                              for( i=12; i<42; i++ ) //24-53
                                buffer[i+12] = (complex_double) clover_pt[i];
                            } else {
                              for( i=0; i<12; i++ ) //0-23
                                buffer[i+12] = buffer[i] = (complex_double) clover_pt[i];
                              for( i=12; i<42; i++ ) //24-53
                                buffer[i+12] = _COMPLEX_double_ZERO;
                            }              
                            for( i=0; i<12; i++ ) //54-65
                              buffer[i+54] = (complex_double) eps_term_pt[i];
#ifdef HAVE_TM
                            if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
                              for(int i=0; i<12; i++) { //0-23
                                buffer[i] += (complex_double) tm_term_pt[i];
                                buffer[i+12] -= (complex_double) tm_term_pt[i];
                              }
                            tm_term_pt += 12;
#endif
                            eps_term_pt += 12;
                            clover_pt += cs;
                            selfcoupling_LU_doublet_decomposition_PRECISION( clover_oo_inv_pt, buffer );
                            clover_oo_inv_pt += 288;
                          }
                        }
                }
#endif  
}

void block_diag_ee_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)  
  int n1 = s->num_block_even_sites, nv = l->num_lattice_site_var;
  clover_PRECISION( eta, phi, &(s->op), start, start+nv*n1, l, threading ); 

  END_UNTHREADED_FUNCTION(threading)
}
                      
// diagonal blocks applied to the odd sites of a block
void block_diag_oo_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int block_num = start/24/(n1+n2);
    //    config_PRECISION clover = s->op.clover_doublet_oo_inv+n1*288+(start/24)*288;
    config_PRECISION clover = s->op.clover_doublet_oo_inv-(block_num+1)*n1*288;
    LU_multiply_PRECISION( eta, phi, clover, n1*24+start, (n1+n2)*24+start );
  } else {
#endif
    if ( g.csw ) {
      int block_num = start/12/(n1+n2);
#ifndef HAVE_TM
      config_PRECISION clover = s->op.clover_oo_inv-(block_num+1)*n1*42;
      LLH_multiply_PRECISION( eta, phi, clover, n1*12+start, (n1+n2)*12+start );
#else
      config_PRECISION clover = s->op.clover_oo_inv-(block_num+1)*n1*72;
      LU_multiply_PRECISION( eta, phi, clover, n1*12+start, (n1+n2)*12+start );
#endif
    } else {
      vector_PRECISION lphi, leta;
      lphi.vector_buffer = phi->vector_buffer+n1*12+start;
      leta.vector_buffer = eta->vector_buffer+n1*12+start;
      config_PRECISION clover = s->op.clover+n1*12+start;
#ifndef HAVE_TM
      for ( i=0; i<12*n2; i++ )
        leta.vector_buffer[i] = lphi.vector_buffer[i]*(clover[i]);
#else
      config_PRECISION tm_term = s->op.tm_term+n1*12+start;
      for ( i=0; i<12*n2; i++ )
        leta.vector_buffer[i] = lphi.vector_buffer[i]*(clover[i]+tm_term[i]);
#endif
    }
#ifdef HAVE_TM1p1
  }
#endif  
  END_UNTHREADED_FUNCTION(threading)
}

// inverted diagonal blocks applied to the odd sites of a block
void block_diag_oo_inv_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s,
                                  level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)
  
  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {

    int block_num = start/24/(n1+n2);
    config_PRECISION clover = s->op.clover_doublet_oo_inv-(block_num+1)*n1*288;
    LU_perform_fwd_bwd_subs_PRECISION( eta, phi, clover, n1*24+start, (n1+n2)*24+start );
  } else {
#endif

    if ( g.csw ) {
      int block_num = start/12/(n1+n2);
#ifndef HAVE_TM
      config_PRECISION clover = s->op.clover_oo_inv-(block_num+1)*n1*42;
      LLH_perform_fwd_bwd_subs_PRECISION( eta, phi, clover, n1*12+start, (n1+n2)*12+start );
#else
      config_PRECISION clover = s->op.clover_oo_inv-(block_num+1)*n1*72;
      LU_perform_fwd_bwd_subs_PRECISION( eta, phi, clover, n1*12+start, (n1+n2)*12+start );
#endif
    } else {
      config_PRECISION clover = s->op.clover+n1*12+start;
      vector_PRECISION lphi, leta;
      lphi.vector_buffer = phi->vector_buffer+n1*12+start;
      leta.vector_buffer = eta->vector_buffer+n1*12+start;
#ifndef HAVE_TM
      for ( i=0; i<12*n2; i++ )
        leta.vector_buffer[i] = lphi.vector_buffer[i]/(clover[i]);
#else
      config_PRECISION tm_term = s->op.tm_term+n1*12+start;
      for ( i=0; i<12*n2; i++ )
        leta.vector_buffer[i] = lphi.vector_buffer[i]/(clover[i]+tm_term[i]);
#endif
    }
#ifdef HAVE_TM1p1
  }
#endif
  
  END_UNTHREADED_FUNCTION(threading)
}


void block_hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
    int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int a1, a2, n1, n2, *length_even = s->dir_length_even, *length_odd = s->dir_length_odd,
      **index = s->oe_index, *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;
  
  config_PRECISION D = s->op.D + (start/nv)*36;
  int i, j, k, *ind;
  config_PRECISION D_pt; 
  buffer_PRECISION lphi = phi->vector_buffer+start, leta = eta->vector_buffer+start;

#ifdef HAVE_TM1p1  
  if ( g.n_flavours == 2 ) {
    complex_PRECISION buf1[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2 = buf1+12;
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[T];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprp_T_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_T_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -T coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprn_T_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_T_PRECISION( buf2, leta+24*j );
    }
  
    // Z direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Z];
      a2 = n1; n2 = a2 + length_odd[Z];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Z]; n1 = a1 + length_odd[Z];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Z]+length_odd[Z];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[Z];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprp_Z_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Z_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Z coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprn_Z_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Z_PRECISION( buf2, leta+24*j );
    }
  
    // Y direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Y];
      a2 = n1; n2 = a2 + length_odd[Y];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Y]; n1 = a1 + length_odd[Y];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Y]+length_odd[Y];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[Y];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprp_Y_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Y_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Y coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprn_Y_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Y_PRECISION( buf2, leta+24*j );
    }
    
    // X direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[X];
      a2 = n1; n2 = a2 + length_odd[X];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[X]; n1 = a1 + length_odd[X];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[X]+length_odd[X];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[X];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprp_X_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_X_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -X coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprn_X_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_X_PRECISION( buf2, leta+24*j );
    }
  } else {
#endif
    complex_PRECISION buf1[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2 = buf1+6;
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
  // "amount" of a block, +T coupling
  ind = index[T];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
    prp_T_PRECISION( buf1, lphi+12*j );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_T_PRECISION( buf2, leta+12*k );
  }
  // "amount" of a block, -T coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
    prn_T_PRECISION( buf1, lphi+12*k );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_T_PRECISION( buf2, leta+12*j );
  }
  
  // Z direction
  if ( amount == _EVEN_SITES ) {
    a1 = 0; n1 = length_even[Z];
    a2 = n1; n2 = a2 + length_odd[Z];
  } else if ( amount == _ODD_SITES ) {
    a1 = length_even[Z]; n1 = a1 + length_odd[Z];
    a2 = 0; n2 = a1;
  } else {
    a1 = 0; n1 = length_even[Z]+length_odd[Z];
    a2 = 0; n2 = n1;
  }
  // "amount" of a block, +Z coupling
  ind = index[Z];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
    prp_Z_PRECISION( buf1, lphi+12*j );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_Z_PRECISION( buf2, leta+12*k );
  }
  // "amount" of a block, -Z coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
    prn_Z_PRECISION( buf1, lphi+12*k );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_Z_PRECISION( buf2, leta+12*j );
  }
  
  // Y direction
  if ( amount == _EVEN_SITES ) {
    a1 = 0; n1 = length_even[Y];
    a2 = n1; n2 = a2 + length_odd[Y];
  } else if ( amount == _ODD_SITES ) {
    a1 = length_even[Y]; n1 = a1 + length_odd[Y];
    a2 = 0; n2 = a1;
  } else {
    a1 = 0; n1 = length_even[Y]+length_odd[Y];
    a2 = 0; n2 = n1;
  }
  // "amount" of a block, +Y coupling
  ind = index[Y];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
    prp_Y_PRECISION( buf1, lphi+12*j );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_Y_PRECISION( buf2, leta+12*k );
  }
  // "amount" of a block, -Y coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
    prn_Y_PRECISION( buf1, lphi+12*k );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_Y_PRECISION( buf2, leta+12*j );
  }
  
  // X direction
  if ( amount == _EVEN_SITES ) {
    a1 = 0; n1 = length_even[X];
    a2 = n1; n2 = a2 + length_odd[X];
  } else if ( amount == _ODD_SITES ) {
    a1 = length_even[X]; n1 = a1 + length_odd[X];
    a2 = 0; n2 = a1;
  } else {
    a1 = 0; n1 = length_even[X]+length_odd[X];
    a2 = 0; n2 = n1;
  }
  // "amount" of a block, +X coupling
  ind = index[X];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
    prp_X_PRECISION( buf1, lphi+12*j );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_X_PRECISION( buf2, leta+12*k );
  }
  // "amount" of a block, -X coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
    prn_X_PRECISION( buf1, lphi+12*k );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_X_PRECISION( buf2, leta+12*j );
  }
#ifdef HAVE_TM1p1
  }
#endif
  END_UNTHREADED_FUNCTION(threading)
}


void block_n_hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
    int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int a1, a2, n1, n2, *length_even = s->dir_length_even, *length_odd = s->dir_length_odd,
      **index = s->oe_index, *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;

  int i, j, k, *ind;
  buffer_PRECISION lphi = phi->vector_buffer+start, leta = eta->vector_buffer+start;
  config_PRECISION D_pt, D = s->op.D + (start/nv)*36;

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    complex_PRECISION buf1[24], *buf2 = buf1+12;
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[T];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprp_T_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_T_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -T coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprn_T_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_T_PRECISION( buf2, leta+24*j );
    }
    
    // Z direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Z];
      a2 = n1; n2 = a2 + length_odd[Z];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Z]; n1 = a1 + length_odd[Z];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Z]+length_odd[Z];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[Z];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprp_Z_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Z_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Z coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprn_Z_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Z_PRECISION( buf2, leta+24*j );
    }
    
    // Y direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Y];
      a2 = n1; n2 = a2 + length_odd[Y];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Y]; n1 = a1 + length_odd[Y];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Y]+length_odd[Y];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[Y];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprp_Y_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Y_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Y coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprn_Y_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Y_PRECISION( buf2, leta+24*j );
    }
    
    // X direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[X];
      a2 = n1; n2 = a2 + length_odd[X];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[X]; n1 = a1 + length_odd[X];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[X]+length_odd[X];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[X];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprp_X_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_X_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -X coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprn_X_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_X_PRECISION( buf2, leta+24*j );
    } 
  } else {
#endif
    complex_PRECISION buf1[12], *buf2 = buf1+6;
    
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[T];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      prp_T_PRECISION( buf1, lphi+12*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      pbp_su3_T_PRECISION( buf2, leta+12*k );
    }
    // "amount" of a block, -T coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      prn_T_PRECISION( buf1, lphi+12*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      pbn_su3_T_PRECISION( buf2, leta+12*j );
    }
    
    // Z direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Z];
      a2 = n1; n2 = a2 + length_odd[Z];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Z]; n1 = a1 + length_odd[Z];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Z]+length_odd[Z];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[Z];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      prp_Z_PRECISION( buf1, lphi+12*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      pbp_su3_Z_PRECISION( buf2, leta+12*k );
    }
    // "amount" of a block, -Z coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      prn_Z_PRECISION( buf1, lphi+12*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      pbn_su3_Z_PRECISION( buf2, leta+12*j );
    }
    
    // Y direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Y];
      a2 = n1; n2 = a2 + length_odd[Y];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Y]; n1 = a1 + length_odd[Y];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Y]+length_odd[Y];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[Y];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      prp_Y_PRECISION( buf1, lphi+12*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      pbp_su3_Y_PRECISION( buf2, leta+12*k );
    }
    // "amount" of a block, -Y coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      prn_Y_PRECISION( buf1, lphi+12*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      pbn_su3_Y_PRECISION( buf2, leta+12*j );
    }
    
    // X direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[X];
      a2 = n1; n2 = a2 + length_odd[X];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[X]; n1 = a1 + length_odd[X];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[X]+length_odd[X];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[X];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      prp_X_PRECISION( buf1, lphi+12*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      pbp_su3_X_PRECISION( buf2, leta+12*k );
    }
    // "amount" of a block, -X coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      prn_X_PRECISION( buf1, lphi+12*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      pbn_su3_X_PRECISION( buf2, leta+12*j );
    }
#ifdef HAVE_TM1p1
  }
#endif
  END_UNTHREADED_FUNCTION(threading)
}


void apply_block_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, int start,
    schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  vector_PRECISION *tmp = s->oe_buf;
  
  block_diag_ee_PRECISION( out, in, start, s, l, threading );
  START_LOCKED_MASTER(threading)
  vector_PRECISION_define( &tmp[0], 0, start + l->num_lattice_site_var*s->num_block_even_sites, start + s->block_vector_size, l );
  END_LOCKED_MASTER(threading)
  block_hopping_term_PRECISION( &tmp[0], in, start, _ODD_SITES, s, l, threading );
  block_diag_oo_inv_PRECISION( &tmp[1], &tmp[0], start, s, l, threading );
  block_n_hopping_term_PRECISION( out, &tmp[1], start, _EVEN_SITES, s, l, threading );
}


void block_solve_oddeven_PRECISION( vector_PRECISION *phi, vector_PRECISION *r, vector_PRECISION *latest_iter,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  vector_PRECISION *tmp = s->oe_buf;
  int end = start+s->block_vector_size;
  
  // odd to even
  vector_PRECISION_copy( &tmp[3], r, start, end, l );
  block_diag_oo_inv_PRECISION( &tmp[2], &tmp[3], start, s, l, no_threading );
  block_n_hopping_term_PRECISION( &tmp[3], &tmp[2], start, _EVEN_SITES, s, l, no_threading );
  local_minres_PRECISION( NULL, &tmp[3], &tmp[2], start, s, l, no_threading );  
  // even to odd
  block_n_hopping_term_PRECISION( &tmp[3], &tmp[2], start, _ODD_SITES, s, l, no_threading );
  block_diag_oo_inv_PRECISION( &tmp[2], &tmp[3], start, s, l, no_threading );
  
  // update phi, latest_iter
  vector_PRECISION_copy( latest_iter, &tmp[2], start, end, l );
  vector_PRECISION_plus( phi, phi, &tmp[2], start, end, l );
  // update r
  vector_PRECISION_copy( r, &tmp[3], start, start+l->num_lattice_site_var*s->num_block_even_sites, l );
  vector_PRECISION_define( r, 0, start+l->num_lattice_site_var*s->num_block_even_sites, end, l );

  END_UNTHREADED_FUNCTION(threading)
}


void block_oddeven_PRECISION_test( level_struct *l, struct Thread *threading ) {
  START_UNTHREADED_FUNCTION(threading)

  schwarz_PRECISION_struct *s = &(l->s_PRECISION);
  
  vector_PRECISION b1, b2, b3, b4, b5;
  PRECISION diff;
  
  vector_PRECISION_init(&b1);
  vector_PRECISION_init(&b2);
  vector_PRECISION_init(&b3);
  vector_PRECISION_init(&b4);
  vector_PRECISION_init(&b5);  

  int vs = s->block_vector_size * s->num_blocks;

  MALLOC( b1.vector_buffer, complex_PRECISION, vs );
  MALLOC( b2.vector_buffer, complex_PRECISION, vs );
  MALLOC( b3.vector_buffer, complex_PRECISION, vs );
  MALLOC( b4.vector_buffer, complex_PRECISION, vs );
  MALLOC( b5.vector_buffer, complex_PRECISION, vs );
  
  vector_PRECISION_define_random( &b1, 0, vs, l );

  for (int i = 0; i< s->num_blocks; i++ ) {
    block_diag_ee_PRECISION( &b2, &b1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
    block_diag_oo_PRECISION( &b2, &b1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
    block_hopping_term_PRECISION( &b2, &b1, s->block[i].start*l->num_lattice_site_var, _FULL_SYSTEM, s, l, no_threading );
    
    block_d_plus_clover_PRECISION( &b3, &b1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
  }
  
  vector_PRECISION_minus( &b3, &b3, &b2, 0, vs, l );
  diff = global_norm_PRECISION( &b3, 0, vs, l, no_threading ) / global_norm_PRECISION( &b2, 0, vs, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of block odd even layout: %le\n", l->depth, diff );
  
  vector_PRECISION_copy( &b4, &b1, 0, s->block_vector_size, l );
  vector_PRECISION_define( &b3, 0, l->num_lattice_site_var*s->num_block_even_sites, s->block_vector_size, l );
  
  block_hopping_term_PRECISION( &b3, &b4, 0, _ODD_SITES, s, l, no_threading );
  block_diag_oo_inv_PRECISION( &b5, &b3, 0, s, l, no_threading );
  vector_PRECISION_plus( &b4, &b4, &b5, l->num_lattice_site_var*s->num_block_even_sites, s->block_vector_size, l );
  
  apply_block_schur_complement_PRECISION( &b3, &b4, 0, s, l, no_threading );
  block_diag_oo_PRECISION( &b3, &b4, 0, s, l, no_threading );
  
  block_diag_oo_inv_PRECISION( &b5, &b3, 0, s, l, no_threading );
  block_hopping_term_PRECISION( &b3, &b5, 0, _EVEN_SITES, s, l, no_threading );
  
  vector_PRECISION_minus( &b3, &b2, &b3, 0, s->block_vector_size, l );
  diff = global_norm_PRECISION( &b3, 0, s->block_vector_size, l, no_threading ) / global_norm_PRECISION( &b2, 0, s->block_vector_size, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of block odd even schur complement: %le\n", l->depth, diff );
  
  FREE( b1.vector_buffer, complex_PRECISION, vs );
  FREE( b2.vector_buffer, complex_PRECISION, vs );
  FREE( b3.vector_buffer, complex_PRECISION, vs );
  FREE( b4.vector_buffer, complex_PRECISION, vs );
  FREE( b5.vector_buffer, complex_PRECISION, vs );

  END_UNTHREADED_FUNCTION(threading)
}

void oddeven_PRECISION_test( level_struct *l ) {

/*********************************************************************************
* In the first step this function checks the correctness of the odd even layout. 
* This is done by:
* - Applying D_W in odd-even structure to a vector.
* - Applying D_W directly to the same vector.
* - Compare solutions ( Difference should be close to 0 ).
*********************************************************************************/  
  
  vector_double d[3];
  vector_PRECISION f[5];
  double diff;
  
  for(int i=0; i<3; i++){
    vector_double_init( &d[i] );
    vector_double_alloc( &d[i], _INNER, 1, l, no_threading );
  }

  for(int i=0; i<5; i++){                                                                 
    vector_PRECISION_init( &f[i] );                                                          
    vector_PRECISION_alloc( &f[i], _INNER, 1, l, no_threading );                             
  } 

  vector_double_define_random( &d[0], 0, l->inner_vector_size, l ); 
  serial_to_oddeven_PRECISION( &f[0], &d[0], l, no_threading );
   
  diag_ee_PRECISION( &f[1], &f[0], &(l->oe_op_PRECISION), l, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var );
  diag_oo_PRECISION( &f[1], &f[0], &(l->oe_op_PRECISION), l, no_threading );
  
  hopping_term_PRECISION( &f[1], &f[0], &(l->oe_op_PRECISION), _FULL_SYSTEM, l, no_threading );
  
  d_plus_clover_double( &d[1], &d[0], &(g.op_double), l, no_threading );
  oddeven_to_serial_PRECISION( &d[0], &f[1], l, no_threading );
  
  vector_double_minus( &d[2], &d[0], &d[1], 0, l->num_inner_lattice_sites, l );
  diff = global_norm_double( &d[2], 0, l->num_inner_lattice_sites, l, no_threading )/global_norm_double( &d[0], 0, l->num_inner_lattice_sites, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of odd even layout: %le\n", l->depth, diff );
    
  // --------------
  
  vector_PRECISION_copy( &f[3], &f[0], 0, l->inner_vector_size, l );
  diag_oo_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), l, no_threading );
  diag_oo_inv_PRECISION( &f[3], &f[2], &(l->oe_op_PRECISION), l, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size );
  vector_PRECISION_minus( &f[3], &f[3], &f[0], 0, l->inner_vector_size, l );

  diff = (PRECISION) (global_norm_PRECISION( &f[3], 0, l->inner_vector_size, l, no_threading )/global_norm_PRECISION( &f[0], 0, l->inner_vector_size, l, no_threading ));
  
  test0_PRECISION("depth: %d, correctness of odd even diagonal term: %le\n", l->depth, diff );
    
  // transformation part
  vector_PRECISION_copy( &f[3], &f[0], 0, l->inner_vector_size, l );
  // even to odd
  // set odd part of f3 to 0. 
  vector_PRECISION_define( &f[2], 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
  
  hopping_term_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), _ODD_SITES, l, no_threading );
  diag_oo_inv_PRECISION( &f[4], &f[2], &(l->oe_op_PRECISION), l, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size );
  vector_PRECISION_plus( &f[3], &f[3], &f[4], l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
  
  // block diagonal part
  apply_schur_complement_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), l, no_threading );
  diag_oo_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), l, no_threading );
  // back transformation part
  diag_oo_inv_PRECISION( &f[4], &f[3], &(l->oe_op_PRECISION), l, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size );
  hopping_term_PRECISION( &f[2], &f[4], &(l->oe_op_PRECISION), _EVEN_SITES, l, no_threading );
  
  vector_PRECISION_minus( &f[0], &f[1], &f[2], 0, l->inner_vector_size, l );
  diff = (PRECISION) (global_norm_PRECISION( &f[0], 0, l->inner_vector_size, l, no_threading )/global_norm_PRECISION( &f[1], 0, l->inner_vector_size, l, no_threading ));
  
  test0_PRECISION("depth: %d, correctness of odd even schur complement: %le\n", l->depth, diff );
  
  for(int i=0; i<3; i++)
    vector_double_free( &d[i], l, no_threading );

  for(int i=0; i<5; i++)
    vector_PRECISION_free( &f[i], l, no_threading );
}
