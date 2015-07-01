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

#include "BlackScholes.h"
#include "complex.h"

namespace Pricer {

	CBlackScholes::CBlackScholes(void) { }
	CBlackScholes::~CBlackScholes(void) { }

	// Waissi and Rossin normal cdf approximation
	double CBlackScholes::normcdf(double z) {
		double pi = 3.141592653589793;
		double b1 = -0.0004406;
		double b2 =  0.0418198;
		double b3 =  0.9;
		return 1.0 / (1.0 + exp(-sqrt(pi)*(b1*pow(z,5.0) + b2*pow(z,3.0) + b3*z)));
	}

	// Black Scholes price
	double CBlackScholes::BSPrice(double S,double K,double r,double q,double v,double T,char PutCall) {
		double d1 = (log(S/K) + (r-q+v*v/2.0)*T)/v/sqrt(T);
		double d2 = d1 - v*sqrt(T);
		double BSCall = S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
		if (PutCall=='C')
			return BSCall;
		else
			return BSCall - S*exp(-q*T) + K*exp(-r*T);
	}

	// Black Scholes characteristic function
	utils::complex CBlackScholes::BlackScholesCF(utils::complex phi,double S0,double r,double v,double T) {
		double u  = log(S0) + (r-v*v/2.0)*T;
		double s2 = v*v*T;
		utils::complex i(0.0,1.0);
		return exp(i*u*phi - 0.5*s2*phi*phi);
	}

}