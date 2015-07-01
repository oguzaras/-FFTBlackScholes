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

#define _in const
#define _out 

enum platforme {
	CPU,
	GPU
};

#ifdef __cplusplus
extern "C" {
#endif

	struct OptionPlan;

	void IntegGridGenerate(double* _v, const double& _eta, const int& _N);

	void LogStrikesGridGenerat(double* _v, const double& _lambdainc, const double& _b, const double& _s0, const int& _N);

	void launch_BlackScholesFFT( const int& N,const double& uplimit, const double& S0, const double& r,const double& q,const double& tau,const double& sigma, 
		const double& alpha, OptionPlan& output);

	void launch_InitGraphicCard();

#ifdef __cplusplus
}
#endif