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

#pragma once

#include <math.h>
#include "complex.h"

using namespace std;

namespace Pricer	{

	class CBlackScholes {
	public:
		CBlackScholes(void);
		~CBlackScholes(void);

		// Waissi and Rossin normal cdf approximation
		double normcdf(double z);

		// Black Scholes price
		double BSPrice(double S,double K,double r,double q,double v,double T,char PutCall);

		// Black Scholes characteristic function
		utils::complex BlackScholesCF(utils::complex phi,double S0,double r,double v,double T);

	};

}