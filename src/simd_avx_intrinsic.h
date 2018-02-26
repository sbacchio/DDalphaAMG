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

#ifndef SIMD_AVX_INTRINSIC_HEADER
#define SIMD_AVX_INTRINSIC_HEADER

#include "immintrin.h"
#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"

#define SIMD               _AVX
#define SIMD_LENGTH_float  8
#define SIMD_LENGTH_double 4
#define mm_FOR_float(e)  { e e e e  e e e e }
#define mm_FOR_double(e) { e e e e }

#define mm_float  __m256
#define mm_double __m256d

#define mm_mul_float  _mm256_mul_ps
#define mm_mul_double _mm256_mul_pd
#define mm_add_float  _mm256_add_ps
#define mm_add_double _mm256_add_pd
#define mm_sub_float  _mm256_sub_ps
#define mm_sub_double _mm256_sub_pd
#define mm_and_float  _mm256_and_ps
#define mm_and_double _mm256_and_pd

#define mm_setzero_float   _mm256_setzero_ps
#define mm_setzero_double  _mm256_setzero_pd
#define mm_setr_float      _mm256_setr_ps
#define mm_setr_double     _mm256_setr_pd
#define mm_set1_float      _mm256_set1_ps
#define mm_set1_double     _mm256_set1_pd
#define mm_load_float      _mm256_load_ps
#define mm_load_double     _mm256_load_pd
#define mm_unpacklo_float  _mm256_unpacklo_ps
#define mm_unpacklo_double _mm256_unpacklo_pd
#define mm_unpackhi_float  _mm256_unpackhi_ps
#define mm_unpackhi_double _mm256_unpackhi_pd
#define mm_store_float    _mm256_store_ps
#define mm_store_double   _mm256_store_pd

#ifdef _FMA_

#define mm_fmadd_float   _mm256_fmadd_ps
#define mm_fmadd_double  _mm256_fmadd_pd
#define mm_fnmadd_float  _mm256_fnmadd_ps
#define mm_fnmadd_double _mm256_fnmadd_pd
#define mm_fmsub_float   _mm256_fmsub_ps
#define mm_fmsub_double  _mm256_fmsub_pd
#define mm_fnmsub_float  _mm256_fnmsub_ps
#define mm_fnmsub_double _mm256_fnmsub_pd

#endif

// Load even components
static inline mm_float mm_seti_float( float *data, const int i ) {
  return mm_setr_float( data[0*i], data[1*i], data[2*i], data[3*i], data[4*i], data[5*i], data[6*i], data[7*i] );
}
static inline mm_double mm_seti_double( double *data, const int i ) {
  return mm_setr_double( data[0*i], data[1*i], data[2*i], data[3*i] );
}

// Loading 6 time the same component and then jumping 12 components 
static inline void mm_set1_6times_float( float *data, mm_float *pack1of3, mm_float *pack2of3,
                                          mm_float *pack3of3, const int skip ) {
  *pack1of3 = mm_setr_float( data[0*i+0*skip], data[1*i+0*skip], data[2*i+0*skip], data[3*i+0*skip],
                             data[4*i+0*skip], data[5*i+0*skip], data[0*i+1*skip], data[1*i+1*skip] );
  *pack2of3 = mm_setr_float( data[2*i+1*skip], data[3*i+1*skip], data[4*i+1*skip], data[5*i+1*skip],
                             data[0*i+2*skip], data[1*i+2*skip], data[2*i+2*skip], data[3*i+2*skip] );
  *pack3of3 = mm_setr_float( data[4*i+2*skip], data[5*i+2*skip], data[0*i+3*skip], data[1*i+3*skip],
                             data[2*i+3*skip], data[3*i+3*skip], data[4*i+3*skip], data[5*i+3*skip] );
}
static inline void mm_loadi_6times_double( double *data, mm_double *pack1of3, mm_double *pack2of3,
                                           mm_double *pack3of3, const int i, const int skip ) {
  *pack1of3 = mm_setr_double( data[0*i+0*skip], data[1*i+0*skip], data[2*i+0*skip], data[3*i+0*skip] );
  *pack2of3 = mm_setr_double( data[4*i+0*skip], data[5*i+0*skip], data[0*i+1*skip], data[1*i+1*skip] );
  *pack3of3 = mm_setr_double( data[2*i+1*skip], data[3*i+1*skip], data[4*i+1*skip], data[5*i+1*skip] );
}

static inline mm_float mm_set_from_list_float( float *data, float *alpha, int *list ) {
  return mm_setr_float( alpha[0]*data[list[0]], alpha[1]*data[list[1]], alpha[2]*data[list[2]], alpha[3]*data[list[3]],
                        alpha[4]*data[list[4]], alpha[5]*data[list[5]], alpha[6]*data[list[6]], alpha[7]*data[list[7]] );
}
static inline mm_double mm_set_from_list_double( double *data, double *alpha, int *list ) {
  return mm_setr_double( alpha[0]*data[list[0]], alpha[1]*data[list[1]], alpha[2]*data[list[2]], alpha[3]*data[list[3]] );
}

// Sum all components of mm_PRECISION
static inline float mm_reduce_add_float( mm_float v) {
  __m128 vlow  = _mm256_castps256_ps128(v);
  __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
  vlow         = _mm_add_ps(vlow, vhigh);     // add the low 128
  // same of SSE
  __m128 shuf  = _mm_movehdup_ps(v);          // broadcast elements 3,1 to 2,0
  __m128 sums  = _mm_add_ps(v, shuf);
  shuf         = _mm_movehl_ps(shuf, sums);   // high half -> low half
  sums         = _mm_add_ss(sums, shuf);
  return       _mm_cvtss_f32(sums);
}
static inline double mm_reduce_add_double( mm_double v ) {
  __m128d vlow  = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);
  vlow          = _mm_add_pd(vlow, vhigh);
  // same of SSE
  double tmp;
  _mm_storeh_pd(&tmp, vlow);        // store the high half
  return _mm_cvtsd_f64(vlow) + tmp; // cast the low half and sum
}

// Transpose a block of SIMD_LENGTH * SIMD_LENGTH
static inline void mm_transpose_float( mm_float *data ) {
  mm_float __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
  mm_float __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7;
  __t0 = _mm256_unpacklo_ps(data[0], data[1]);
  __t1 = _mm256_unpackhi_ps(data[0], data[1]);
  __t2 = _mm256_unpacklo_ps(data[2], data[3]);
  __t3 = _mm256_unpackhi_ps(data[2], data[3]);
  __t4 = _mm256_unpacklo_ps(data[4], data[5]);
  __t5 = _mm256_unpackhi_ps(data[4], data[5]);
  __t6 = _mm256_unpacklo_ps(data[6], data[7]);
  __t7 = _mm256_unpackhi_ps(data[6], data[7]);
  __tt0 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(1,0,1,0));
  __tt1 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(3,2,3,2));
  __tt2 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(1,0,1,0));
  __tt3 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(3,2,3,2));
  __tt4 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(1,0,1,0));
  __tt5 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(3,2,3,2));
  __tt6 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(1,0,1,0));
  __tt7 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(3,2,3,2));
  data[0] = _mm256_permute2f128_ps(__tt0, __tt4, 0x20);
  data[1] = _mm256_permute2f128_ps(__tt1, __tt5, 0x20);
  data[2] = _mm256_permute2f128_ps(__tt2, __tt6, 0x20);
  data[3] = _mm256_permute2f128_ps(__tt3, __tt7, 0x20);
  data[4] = _mm256_permute2f128_ps(__tt0, __tt4, 0x31);
  data[5] = _mm256_permute2f128_ps(__tt1, __tt5, 0x31);
  data[6] = _mm256_permute2f128_ps(__tt2, __tt6, 0x31);
  data[7] = _mm256_permute2f128_ps(__tt3, __tt7, 0x31);
}
static inline void mm_transpose_double( mm_double *data)
{
   mm_double tmp[4];

   tmp[0] = _mm256_unpacklo_pd( data[0], data[1] );
   tmp[1] = _mm256_unpacklo_pd( data[2], data[3] );
   tmp[2] = _mm256_unpackhi_pd( data[0], data[1] );
   tmp[3] = _mm256_unpackhi_pd( data[2], data[3] );
   //TODO
   data[0] = _mm256_movelh_pd( tmp[0], tmp[1] );
   data[1] = _mm256_movehl_pd( tmp[1], tmp[0] );
   data[2] = _mm256_movelh_pd( tmp[2], tmp[3] );
   data[3] = _mm256_movehl_pd( tmp[3], tmp[2] );
}

#endif
