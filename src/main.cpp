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


//#include "stdafx.h"
#include "BlackScholes.h"
#include <iostream>
#include <iomanip>

#include "fft_kernel.h"

using namespace std;

// PLATFORM = GPU or CPU
#define PLATFORM GPU

int main() {
	double S0 = 10.0;
	double r = 0.0;
	double q = 0.0;
	double tau = 0.5;
	double sigma = 0.3;
	double alpha = 1.75;
	int Ngrid = pow(2,9.0);
	double uplimit = 200.0;

	OptionPlan output;

	//----- FFT Prices
	//	template argument : CPU or GPU
	Pricer::fft_kernel<PLATFORM>::InitGraphicCard();
	Pricer::fft_kernel<PLATFORM>::exec(Ngrid,uplimit,S0,r,q,tau,sigma,alpha, output);	// FFT Black Scholes prices	

	std::vector<double> K = output.K;
	std::vector<double> CallFFT = output.CallFFT;
	double eta = output.eta;
	double lambda = output.lambdainc;
	int N = K.size();
	//----- Black Scholes prices
	Pricer::CBlackScholes BS;
	std::vector<double> BSCall(N,0.0);
	for (int k=0; k<N; k++)
		BSCall[k] = BS.BSPrice(S0,K[k],r,q,sigma,tau,'C');


	printf("FFT Using %d grid points \n",Ngrid);
	printf("Integration grid increment (eta)    %6.4f \n",eta);
	printf("Log strike grid increment (lambda)  %6.4f \n",lambda);
	printf("-----------------------------\n");
	printf(" Strike   FFTPrice   BSPrice \n");
	printf("-----------------------------\n");
	for (int k=0; k<N; k++) {
		if ((K[k] > 8.0) & (K[k] < 12.0)) 
			printf(" %5.2f    %6.4f     %6.4f \n",K[k],CallFFT[k],BSCall[k]);
	}
	printf("-----------------------------\n");
}

