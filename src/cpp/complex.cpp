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

#include "complex.h"

namespace utils {

	utils::complex exp(const utils::complex& _z)	{
	utils::complex res;
	double t = (double)_z.m_re;
	t		 = expf (t);
	res.m_re = t*cos(_z.m_im);
	res.m_im = t*sin(_z.m_im);
	return res;
}

}