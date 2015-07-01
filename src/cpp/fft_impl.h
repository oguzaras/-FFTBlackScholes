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

#ifndef _FFT_H_
#define _FFT_H_

#include "complex.h"

namespace Pricer {

	class CFFT
	{
	public:
		//   FORWARD FOURIER TRANSFORM
		static bool Forward(const utils::complex *const Input, utils::complex *const Output, const unsigned int N);

		//   FORWARD FOURIER TRANSFORM, INPLACE VERSION
		static bool Forward(utils::complex *const Data, const unsigned int N);

		//   INVERSE FOURIER TRANSFORM
		static bool Inverse(const utils::complex *const Input, utils::complex *const Output, const unsigned int N, const bool Scale = true);

		//   INVERSE FOURIER TRANSFORM, INPLACE VERSION
		static bool Inverse(utils::complex *const Data, const unsigned int N, const bool Scale = true);

	protected:
		//   Addapt input data to the algorithm
		static void AddaptData(const utils::complex *const Input, utils::complex *const Output, const unsigned int N);
		static void AddaptData(utils::complex *const Data, const unsigned int N);

		//   FFT implementation
		static void Transform(utils::complex *const Data, const unsigned int N, const bool Inverse = false);

		//   Scaling of inverse FFT result
		static void Scale(utils::complex *const Data, const unsigned int N);
	};
}
#endif
