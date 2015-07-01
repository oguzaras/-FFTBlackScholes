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

#include "kernels_wrapper.h"
#include <thrust\device_vector.h>
#include <thrust\host_vector.h>
#include <cufft.h>

#include "../common/common.h"
#include "../../src/cuda/cuDblComplex.cuh"

#define _in const
#define _out 

#define BLK_SIZE 256


__device__ cuDoubleComplex BlackScholesCF(cuDoubleComplex phi,double S0,double r,double v,double T) {
	double u  = log(S0) + (r-v*v/2.0)*T;
	double s2 = v*v*T;
	cuDoubleComplex i = make_cuDoubleComplex(0.0,1.0);
	return exp(i*u*phi - 0.5*s2*phi*phi);
}

__global__ void compute_input_data(cuDoubleComplex* x, int N,double lambdainc, double eta, double S0, double r, double q, double tau, double sigma, double alpha, double s0 ){

	double b = double(N)*lambdainc/2.0;
	cuDoubleComplex i= make_cuDoubleComplex(0.0,1.0);
	cuDoubleComplex psi,phi;

	// to optimize
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	bool is_bound = (idx == 0 || idx == N-1);
	if( idx < N ){
		double v = idx * eta;
		psi = BlackScholesCF(v-(alpha+1.0)*i,S0,r,sigma,tau);
		phi = exp(-r*tau)*psi/(alpha*alpha + alpha - v*v+ i*v*(2.0*alpha+1.0));
		x[idx] = exp(i*(b-s0)*v) * phi * ( is_bound ? 0.5 : 1.0);
	}
}


__global__ void compute_strike(double* out, double _b, double _lambdainc, double _s0, int N){
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if( idx < N )
		out[idx]= exp(-_b + _lambdainc*idx + _s0);
}

void launch_BlackScholesFFT( _in int& N,_in double& uplimit,_in double& S0,_in double& r,_in double& q,_in double& tau,_in double& sigma, 
							_in double& alpha, OptionPlan& output)
{
	// log spot price
	double s0 = log(S0);

	// Specify the increments
	double pi = 3.141592653589793;
	double eta = uplimit/double(N);
	double  lambdainc = 2.0*pi/double(N)/eta;
	
	// Specify the b parameter
	double b = double(N)*lambdainc/2.0;
	
	// Create the strikes
	thrust::device_vector<double> dK(N,1.0);
	dim3 grid_dim((N + BLK_SIZE - 1) / BLK_SIZE); 
	compute_strike<<<grid_dim, BLK_SIZE>>>((double*)thrust::raw_pointer_cast(dK.data()), b, lambdainc, s0, N);

	// Initialize the price vectors;
	thrust::device_vector<cuDoubleComplex> dx(N);
	compute_input_data<<<grid_dim, BLK_SIZE>>>( (cuDoubleComplex*)thrust::raw_pointer_cast(dx.data()),
												 N, lambdainc, eta, S0, r, q, tau, sigma, alpha, s0); 
	// apply fft
	cufftHandle plan;
	cufftResult_t ret =  cufftPlan1d(&plan, N, CUFFT_Z2Z, 10);
	ret = cufftExecZ2Z(plan, (cuDoubleComplex*)thrust::raw_pointer_cast(dx.data()), (cuDoubleComplex*)thrust::raw_pointer_cast(dx.data()), CUFFT_FORWARD);

	thrust::host_vector<cuDoubleComplex> CallFFT(N);
	cudaMemcpy(&CallFFT[0], (cuDoubleComplex*)thrust::raw_pointer_cast(dx.data()), sizeof(cuDoubleComplex)*N, cudaMemcpyDeviceToHost);

	output.CallFFT.resize(N);
	output.K.resize(N);


	// utiliser trust poue -- une foi nettoyage effectuer
	for (int u=0; u<=N-1; u++)
		output.CallFFT[u] = eta*exp(-alpha*(-b + lambdainc*u + s0))/pi * CallFFT[u].x;
	
	cudaMemcpy(&(output.K[0]), (double*)thrust::raw_pointer_cast(dK.data()), sizeof(double)*N, cudaMemcpyDeviceToHost);
	output.lambdainc = lambdainc;
	output.eta = eta;

}

__global__ void vec_add(double* v1, double* v2, double* v3, int N)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if( idx < N )
		v1[idx] = v2[idx] + v3[idx];
}

void launch_InitGraphicCard()
{
	int N = 10e6;
	thrust::device_vector<double> dv(N,3.6);
	thrust::device_vector<double> dw(N,15.3);
	thrust::device_vector<double> dr(N,0.0);
	double* hr = new double[N];
	dim3 grid_dim((N + BLK_SIZE - 1) / BLK_SIZE); 
	vec_add<<<grid_dim, BLK_SIZE>>>((double*)thrust::raw_pointer_cast(dv.data()), (double*)thrust::raw_pointer_cast(dw.data()), (double*)thrust::raw_pointer_cast(dr.data()), N);
	cudaMemcpy(hr, (double*)thrust::raw_pointer_cast(dr.data()), sizeof(double)*N, cudaMemcpyDeviceToHost);
	delete[] hr;
}

// cuDoubleComplex : overload operator definition

__host__ __device__ cuDoubleComplex operator+(cuDoubleComplex _lhs, cuDoubleComplex _rhs)	{
	return cuCadd(_lhs, _rhs);
}
__host__ __device__ cuDoubleComplex operator-(cuDoubleComplex _lhs, cuDoubleComplex _rhs)	{
	return cuCsub(_lhs, _rhs);
}
__host__ __device__ cuDoubleComplex operator*(cuDoubleComplex _lhs, cuDoubleComplex _rhs)	{
	return cuCmul(_lhs, _rhs);
}
__host__ __device__ cuDoubleComplex operator/(cuDoubleComplex _lhs, cuDoubleComplex _rhs)	{
	return cuCdiv(_lhs, _rhs);
}

__host__ __device__ cuDoubleComplex operator+(cuDoubleComplex _lhs, double _rhs)	{
	return make_cuDoubleComplex(_lhs.x + _rhs, _lhs.y);
}
__host__ __device__ cuDoubleComplex operator-(cuDoubleComplex _lhs, double _rhs)	{
	return make_cuDoubleComplex(_lhs.x - _rhs, _lhs.y);
}
__host__ __device__ cuDoubleComplex operator*(cuDoubleComplex _lhs, double _rhs)	{
	return make_cuDoubleComplex(_lhs.x * _rhs, _lhs.y * _rhs);
}
__host__ __device__ cuDoubleComplex operator/(cuDoubleComplex _lhs, double _rhs)	{
	return make_cuDoubleComplex(_lhs.x / _rhs, _lhs.y / _rhs);
}

__host__ __device__ cuDoubleComplex operator+(double _lhs, cuDoubleComplex _rhs)	{
	return _rhs + _lhs;
}
__host__ __device__ cuDoubleComplex operator-(double _lhs, cuDoubleComplex _rhs)	{
	return _lhs + (-1.0*_rhs);
}
__host__ __device__ cuDoubleComplex operator*(double _lhs, cuDoubleComplex _rhs)	{
	return _rhs * _lhs;
}
__host__ __device__ cuDoubleComplex operator/(double _lhs, cuDoubleComplex _rhs)	{
	return cuCdiv(make_cuDoubleComplex(_lhs,0) , _rhs);
}

__host__ __device__ cuDoubleComplex exp(cuDoubleComplex _z)	{
	cuDoubleComplex res;
	double t = exp (cuCreal(_z));
	res.x = t*cos(_z.y);
	res.y = t*sin(_z.y);
	return res;
}