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

#include "fft_impl.h"
#include <math.h>

namespace Pricer {

	bool CFFT::Forward(const utils::complex *const Input, utils::complex *const Output, const unsigned int N)
	{
		if (!Input || !Output || N < 1 || N & (N - 1))
			return false;
		AddaptData(Input, Output, N);
		Transform(Output, N);
		return true;
	}

	bool CFFT::Forward(utils::complex *const Data, const unsigned int N)
	{
		if (!Data || N < 1 || N & (N - 1))
			return false;
		AddaptData(Data, N);
		Transform(Data, N);
		return true;
	}

	bool CFFT::Inverse(const utils::complex *const Input, utils::complex *const Output, const unsigned int N, const bool Scale /* = true */)
	{
		if (!Input || !Output || N < 1 || N & (N - 1))
			return false;
		AddaptData(Input, Output, N);
		Transform(Output, N, true);
		if (Scale)
			CFFT::Scale(Output, N);
		return true;
	}

	bool CFFT::Inverse(utils::complex *const Data, const unsigned int N, const bool Scale /* = true */)
	{
		if (!Data || N < 1 || N & (N - 1))
			return false;
		AddaptData(Data, N);
		Transform(Data, N, true);
		if (Scale)
			CFFT::Scale(Data, N);
		return true;
	}

	//   Addapt input data to the algorithm
	void CFFT::AddaptData(const utils::complex *const Input, utils::complex *const Output, const unsigned int N)
	{
		unsigned int Target = 0;
		for (unsigned int Position = 0; Position < N; ++Position)
		{
			Output[Target] = Input[Position];
			unsigned int Mask = N;
			while (Target & (Mask >>= 1))
					Target &= ~Mask;
			Target |= Mask;
		}
	}

	//   Inplace version 
	void CFFT::AddaptData(utils::complex *const Data, const unsigned int N)
	{
		unsigned int Target = 0;
		for (unsigned int Position = 0; Position < N; ++Position)
		{
			if (Target > Position)
			{
				const utils::complex Temp(Data[Target]);
				Data[Target] = Data[Position];
				Data[Position] = Temp;
			}
			unsigned int Mask = N;
			while (Target & (Mask >>= 1))
					Target &= ~Mask;
			Target |= Mask;
		}
	}

	//   FFT implementation
	void CFFT::Transform(utils::complex *const Data, const unsigned int N, const bool Inverse /* = false */)
	{
		const double pi = Inverse ? 3.14159265358979323846 : -3.14159265358979323846;

		for (unsigned int Step = 1; Step < N; Step <<= 1)
		{
			//   Jump to the next entry of the same transform factor
			const unsigned int Jump = Step << 1;
			//   Angle increment
			const double delta = pi / double(Step);
			//    sin(delta / 2)
			const double Sine = sin(delta * .5);
			//   Multiplier for trigonometric recurrence
			const utils::complex Multiplier(-2. * Sine * Sine, sin(delta));
			//   Start value for transform factor, fi = 0
			utils::complex Factor(1.);
			//   Iteration through groups of different transform factor
			for (unsigned int Group = 0; Group < Step; ++Group)
			{
				//   Iteration within group 
				for (unsigned int Pair = Group; Pair < N; Pair += Jump)
				{
					//   Match position
					const unsigned int Match = Pair + Step;
					//   Second term of two-point transform
					const utils::complex Product(Factor * Data[Match]);
					//   Transform for fi + pi
					Data[Match] = Data[Pair] - Product;
					//   Transform for fi
					Data[Pair] += Product;
				}
				//   Successive transform factor via trigonometric recurrence
				Factor = Multiplier * Factor + Factor;
			}
		}
	}

	//   Scaling of inverse FFT result
	void CFFT::Scale(utils::complex *const Data, const unsigned int N)
	{
		const double Factor = 1. / double(N);
		for (unsigned int Position = 0; Position < N; ++Position)
			Data[Position] *= Factor;
	}
}