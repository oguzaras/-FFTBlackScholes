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

#include "FFTOptionPricing.h"
#include "BlackScholes.h"

#include "fft_impl.h"
#include "complex.h"

namespace Pricer {

	CFastFourierTransform::CFastFourierTransform(void) { }
	CFastFourierTransform::~CFastFourierTransform(void) { }




	void price(std::vector<utils::complex >& x, const std::vector<double>& w, const std::vector<double>& v, int N, Pricer::CBlackScholes& BS,double alpha, double S0, double r, double sigma, double tau, double s0, double b){
		utils::complex i(0.0,1.0);
		utils::complex psi,phi;
		for (int j=0; j<=N-1; j++) {
			psi = BS.BlackScholesCF(v[j]-(alpha+1.0)*i,S0,r,sigma,tau);
			phi = exp(-r*tau)*psi/(alpha*alpha + alpha - v[j]*v[j] + i*v[j]*(2.0*alpha+1.0));
			x[j] = exp(i*(b-s0)*v[j])*phi*w[j];
		}
	}

	void CFastFourierTransform::BlackScholesFFT(int N,double uplimit,double S0,double r,double q,double tau,double sigma,double alpha, OptionPlan& output)
	{
		// log spot price
		double s0 = log(S0);

		// Specify the increments
		double eta = uplimit/double(N);
		double pi = 3.141592653589793;
		double lambdainc = 2.0*pi/double(N)/eta;

		// Initialize and specify the weights for the trapezoidal rule
		std::vector<double> w(N,1.0);
		w[0] = 0.5;
		w[N-1] = 0.5;

		// Specify the b parameter
		double b = double(N)*lambdainc/2.0;

		// Create the grid for the integration
		std::vector<double> v(N,0.0);
		for (int n=0; n<=N-1; n++)
			v[n] = eta*n;

		// Create the grid for the log-strikes
		std::vector<double> k(N,0.0);
		for (int n=0; n<=N-1; n++)
			k[n] = -b + lambdainc*n + s0;

		// Create the strikes
		std::vector<double> K(N,0.0);
		for (int n=0; n<=N-1; n++)
			K[n] = exp(k[n]);

		// Initialize the price vectors;
		std::vector<double> CallFFT(N,0.0);

		// Implement the FFT
		Pricer::CBlackScholes BS;

		utils::complex i(0.0,1.0);
		utils::complex psi,phi;

		std::vector< utils::complex > x(N);
		price(x, w, v, N, BS, alpha, S0, r, sigma, tau, s0, b);

		Pricer::CFFT::Forward(&x[0],N);

		for (int u=0; u<=N-1; u++)
			CallFFT[u] = eta*exp(-alpha*k[u])/pi * x[u].m_re;

		output.CallFFT = CallFFT;
		output.K = K;
		output.lambdainc = lambdainc;
		output.eta = eta;
	}
}