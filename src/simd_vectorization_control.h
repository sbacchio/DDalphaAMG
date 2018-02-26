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

#ifndef SIMD_VECTORIZATION_CONTROL_HEADER
#define SIMD_VECTORIZATION_CONTROL_HEADER

#ifdef        NOT_YET//__AVX__
#include "simd_avx_intrinsic.h"
#elif defined SSE //__SSE__
#include "simd_sse_intrinsic.h"
#endif

#ifdef        SIMD

#define OPTIMIZED_COARSE_NEIGHBOR_COUPLING_float
#define OPTIMIZED_COARSE_SELF_COUPLING_float
#define OPTIMIZED_INTERPOLATION_OPERATOR_float
#define OPTIMIZED_INTERPOLATION_SETUP_float
#define OPTIMIZED_NEIGHBOR_COUPLING_double
#define OPTIMIZED_NEIGHBOR_COUPLING_float
//#define OPTIMIZED_SELF_COUPLING_double
#define OPTIMIZED_SELF_COUPLING_float
#define OPTIMIZED_LINALG_float
#define OPTIMIZED_LINALG_double

#define OPERATOR_COMPONENT_OFFSET_float  (SIMD_LENGTH_float *((l->num_eig_vect+SIMD_LENGTH_float -1)/SIMD_LENGTH_float ))
#define OPERATOR_COMPONENT_OFFSET_double (SIMD_LENGTH_double*((l->num_eig_vect+SIMD_LENGTH_double-1)/SIMD_LENGTH_double))

#define OPERATOR_TYPE_float float
#define OPERATOR_TYPE_double double

#endif

#ifndef       __FMA__
// a*b + c
static inline mm_double mm_fmadd_double( mm_double a, mm_double b, mm_double c ) {
  return mm_add_double( mm_mul_double( a, b ), c );
}
static inline mm_float mm_fmadd_float( mm_float a, mm_float b, mm_float c ) {
  return mm_add_float( mm_mul_float( a, b ), c );
}

// -a*b + c
static inline mm_double mm_fnmadd_double( mm_double a, mm_double b, mm_double c ) {
  return mm_sub_double( c, mm_mul_double( a, b ) );
}
static inline mm_float mm_fnmadd_float( mm_float a, mm_float b, mm_float c ) {
  return mm_sub_float( c, mm_mul_float( a, b ) );
}

// a*b - c
static inline mm_double mm_fmsub_double( mm_double a, mm_double b, mm_double c ) {
  return mm_sub_double( mm_mul_double( a, b ), c );
}
static inline mm_float mm_fmsub_float( mm_float a, mm_float b, mm_float c ) {
  return mm_sub_float( mm_mul_float( a, b ), c );
}

// res = -a*b - c
static inline mm_double mm_fnmsub_double( mm_double a, mm_double b, mm_double c ) {
  mm_double na = mm_sub_double( mm_setzero_double(), a );
  return mm_sub_double( mm_mul_double( na, b ), c );
}
static inline mm_float mm_fnmsub_float( mm_float a, mm_float b, mm_float c ) {
  mm_float na = mm_sub_float( mm_setzero_float(), a );
  return mm_sub_float( mm_mul_float( na, b ), c );
}
#endif

#endif // SIMD_VECTORIZATION_CONTROL_HEADER
