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
 *	  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <stdlib.h>
#include <vector>
#include "BlackScholes.h"

#include "kernels_wrapper.h"
#include "FFTOptionPricing.h"

namespace Pricer {

	template< int T>
	struct fft_kernel {
		static void exec(_in int& N,_in double& uplimit,_in double& S0,_in double& r,_in double& q,_in double& tau,_in double& sigma, 
			_in double& alpha, _out OptionPlan& output)
		{
		}
		static void InitGraphicCard(){	}

	};


	template< >
	struct fft_kernel<CPU> {
		static void exec(_in int& N,_in double& uplimit,_in double& S0,_in double& r,_in double& q,_in double& tau,_in double& sigma, 
			_in double& alpha, _out OptionPlan& output)
		{
			CFastFourierTransform BFFT;
			BFFT.BlackScholesFFT(N, uplimit, S0, r, q, tau, sigma, alpha, output);
		}
		static void InitGraphicCard(){	}
	};



	template< >
	struct fft_kernel<GPU> {
		static void exec(_in int& N,_in double& uplimit,_in double& S0,_in double& r,_in double& q,_in double& tau,_in double& sigma, 
			_in double& alpha, _out OptionPlan& output)
		{
			launch_BlackScholesFFT(N, uplimit, S0, r, q, tau, sigma, alpha, output);
		}
		static void InitGraphicCard()
		{
			launch_InitGraphicCard();
		}

	};
}
