/*
	Numerical Integration by Gauss-Legendre Quadrature Formulas of high orders.
	High-precision abscissas and weights are used.

	Project homepage: http://www.holoborodko.com/pavel/?page_id=679
	Contact e-mail:   pavel@holoborodko.com

	Copyright (c)2007-2008 Pavel Holoborodko

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

	Contributors
	Konstantin Holoborodko - Optimization of Legendre polynomial computing.
*/

#ifndef __GAUSS_LEGENDRE_H__
#define __GAUSS_LEGENDRE_H__

#ifdef __cplusplus
extern "C"
{
#endif

	/* Numerical computation of int(f(x),x=a..b) by Gauss-Legendre n-th order high precision quadrature 
		[in]n     - quadrature order
		[in]f     - integrant
		[in]data  - pointer on user-defined data which will 
				   be passed to f every time it called (as second parameter).
		[in][a,b] - interval of integration
	   
		return:
	         -computed integral value or -1.0 if n order quadrature is not supported
	*/

	double gauss_legendre(int n, double (*f)(double,void*), void* data, double a, double b);

	/* Computing of abscissas and weights for Gauss-Legendre quadrature for any(reasonable) order n
		[in] n   - order of quadrature
		[in] eps - required precision (must be eps>=macheps(double), usually eps = 1e-10 is ok)
		[out]x   - abscisass, size = (n+1)>>1
		[out]w   - weights, size = (n+1)>>1 
	*/
	void gauss_legendre_tbl(int n, double* x, double* w, double eps);

#ifdef __cplusplus
}
#endif

#endif /* __GAUSS_LEGENDRE_H__ */


