/*
 * Copyright (C) 2016 Simone Bacchio.
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

#ifndef SIMD_SSE_INTRINSIC_HEADER
#define SIMD_SEE_INTRINSIC_HEADER

#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"

#define SIMD               _SSE
#define SIMD_LENGTH_float  4
#define SIMD_LENGTH_double 2

#define mm_float  __m128
#define mm_double __m128d

#define mm_mul_float  _mm_mul_ps
#define mm_mul_double _mm_mul_pd
#define mm_add_float  _mm_add_ps
#define mm_add_double _mm_add_pd
#define mm_sub_float  _mm_sub_ps
#define mm_sub_double _mm_sub_pd
#define mm_and_float  _mm_and_ps
#define mm_and_double _mm_and_pd

#define mm_setzero_float   _mm_setzero_ps
#define mm_setzero_double  _mm_setzero_pd
#define mm_setr_float      _mm_setr_ps
#define mm_setr_double     _mm_setr_pd
#define mm_set1_float      _mm_set1_ps
#define mm_set1_double     _mm_set1_pd
#define mm_load_float      _mm_load_ps
#define mm_load_double     _mm_load_pd
#define mm_unpacklo_float  _mm_unpacklo_ps
#define mm_unpacklo_double _mm_unpacklo_pd
#define mm_unpackhi_float  _mm_unpackhi_ps
#define mm_unpackhi_double _mm_unpackhi_pd
#define mm_store_float    _mm_store_ps
#define mm_store_double   _mm_store_pd

#ifdef __FMA__

#include "immintrin.h"

#define mm_fmadd_float   _mm_fmadd_ps
#define mm_fmadd_double  _mm_fmadd_pd
#define mm_fnmadd_float  _mm_fnmadd_ps
#define mm_fnmadd_double _mm_fnmadd_pd
#define mm_fmsub_float   _mm_fmsub_ps
#define mm_fmsub_double  _mm_fmsub_pd
#define mm_fnmsub_float  _mm_fnmsub_ps
#define mm_fnmsub_double _mm_fnmsub_pd

#endif

// set components from data, with increment i
static inline mm_float mm_seti_float( float *data, const int i ) {
  return mm_setr_float( data[0*i], data[1*i], data[2*i], data[3*i] );
}
static inline mm_double mm_seti_double( double *data, const int i ) {
  return mm_setr_double( data[0*i], data[1*i] );
}

// Loading 6 times the same component and then jumping *skip* components 
static inline void mm_loadi_6times_float( float *data, mm_float *pack1of3, mm_float *pack2of3,
                                         mm_float *pack3of3, const int i, const int skip ) {
  *pack1of3 = mm_setr_float( data[0*i+0*skip], data[1*i+0*skip], data[2*i+0*skip], data[3*i+0*skip] );
  *pack2of3 = mm_setr_float( data[4*i+0*skip], data[5*i+0*skip], data[0*i+1*skip], data[1*i+1*skip] );
  *pack3of3 = mm_setr_float( data[2*i+1*skip], data[3*i+1*skip], data[4*i+1*skip], data[5*i+1*skip] );
}
static inline void mm_loadi_6times_double( double *data, mm_double *pack1of3, mm_double *pack2of3,
                                           mm_double *pack3of3, const int i, const int skip ) {
  *pack1of3 = mm_setr_double( data[0*i+0*skip], data[1*i+0*skip] );
  *pack2of3 = mm_setr_double( data[2*i+0*skip], data[3*i+0*skip] );
  *pack3of3 = mm_setr_double( data[4*i+0*skip], data[5*i+0*skip] );
}

// Set from list and scaling with alpha 
static inline mm_float mm_set_from_list_float( float *data, float *alpha, int *list ) {
  return mm_setr_float( alpha[0]*data[list[0]], alpha[1]*data[list[1]], alpha[2]*data[list[2]], alpha[3]*data[list[3]] );
}
static inline mm_double mm_set_from_list_double( double *data, double *alpha, int *list ) {
  return mm_setr_double( alpha[0]*data[list[0]], alpha[1]*data[list[1]] );
}

// Sum all components of mm_PRECISION
static inline float mm_reduce_add_float( mm_float v ) {
  mm_float shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
  mm_float sums = _mm_add_ps(v, shuf);
  shuf          = _mm_movehl_ps(shuf, sums); // high half -> low half
  sums          = _mm_add_ss(sums, shuf);
  return        _mm_cvtss_f32(sums);
}
static inline double mm_reduce_add_double( mm_double v ) {
  double tmp;
  _mm_storeh_pd(&tmp, v);        // store the high half
  return _mm_cvtsd_f64(v) + tmp; // cast the low half and sum
}

// Transpose a block of SIMD_LENGTH * SIMD_LENGTH
static inline void mm_transpose_float( mm_float *data ) {
  _MM_TRANSPOSE4_PS(data[0],data[1],data[2],data[3]);
}
static inline void mm_transpose_double( mm_double *data ) {
  double tmp01, tmp10 = _mm_cvtsd_f64(data[1]);
  _mm_storeh_pd(&tmp01, data[0]);
  _mm_loadl_pd(data[1], &tmp01);
  _mm_loadh_pd(data[0], &tmp10);
}

#endif 
