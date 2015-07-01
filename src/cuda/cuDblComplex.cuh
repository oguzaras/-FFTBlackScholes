/*
 *  Copyright (C) 2015  Oguz ARAS
 *
 *	This file is part of FFTBlackScholes.
 *
 *    FFTBlackScholes is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    FFTBlackScholes is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CU_COMPLEX_CUSTOM_
#define _CU_COMPLEX_CUSTOM_


#include <cuComplex.h>

// cuDoubleComple: overload operator
__host__ __device__ cuDoubleComplex operator+(cuDoubleComplex _lhs, cuDoubleComplex _rhs);  
__host__ __device__ cuDoubleComplex operator-(cuDoubleComplex _lhs, cuDoubleComplex _rhs);
__host__ __device__ cuDoubleComplex operator*(cuDoubleComplex _lhs, cuDoubleComplex _rhs);
__host__ __device__ cuDoubleComplex operator/(cuDoubleComplex _lhs, cuDoubleComplex _rhs);
__host__ __device__ cuDoubleComplex operator+(cuDoubleComplex _lhs, double _rhs)	;
__host__ __device__ cuDoubleComplex operator-(cuDoubleComplex _lhs, double _rhs)	;
__host__ __device__ cuDoubleComplex operator*(cuDoubleComplex _lhs, double _rhs)	;
__host__ __device__ cuDoubleComplex operator/(cuDoubleComplex _lhs, double _rhs)	;
__host__ __device__ cuDoubleComplex operator+(double _lhs, cuDoubleComplex _rhs)	;
__host__ __device__ cuDoubleComplex operator-(double _lhs, cuDoubleComplex _rhs)	;
__host__ __device__ cuDoubleComplex operator*(double _lhs, cuDoubleComplex _rhs)	;
__host__ __device__ cuDoubleComplex operator/(double _lhs, cuDoubleComplex _rhs)	;


// mathematical functions
__host__ __device__ cuDoubleComplex exp(cuDoubleComplex _z);

#endif //_CU_COMPLEX_CUSTOM_