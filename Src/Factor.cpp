/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

//////////////////////
// Polynomial Roots //
//////////////////////
#include <cmath>

#include "Factor.h"

static double const SQRT_3 = 1.7320508075688772935;

template<>
std::vector<std::complex<double> > Factor<1>(double const cs[2], double EPS) {
	double a0 = cs[0];
	double a1 = cs[1];
	std::vector<std::complex<double> > roots;

	if(fabs(a1) <= EPS) return roots;
	roots.push_back(-a0 / a1);

	return roots;
}

template<>
std::vector<std::complex<double> > Factor<2>(double const cs[3], double EPS) {
	double a0 = cs[0];
	double a1 = cs[1];
	double a2 = cs[2];
	std::vector<std::complex<double> > roots;

	if(fabs(a2) <= EPS) return Factor<1>(cs, EPS);

	double d = a1 * a1 - 4 * a0 * a2;
	a1 /= 2 * a2;
	if(d < 0) {
		d = std::sqrt(-d) / (2 * a2);
		roots.push_back(std::complex<double>(-a1, -d));
		roots.push_back(std::complex<double>(-a1, d));
	} else {
		d = std::sqrt(d) / (2 * a2);
		roots.push_back(-a1 - d);
		roots.push_back(-a1 + d);
	}

	return roots;
}

// Solution taken from: http://mathworld.wolfram.com/CubicFormula.html
// and http://www.csit.fsu.edu/~burkardt/f_src/subpak/subpak.f90
template<>
std::vector<std::complex<double> > Factor<3>(double const cs[4], double EPS) {
	double a0 = cs[0];
	double a1 = cs[1];
	double a2 = cs[2];
	double a3 = cs[3];
	std::vector<std::complex<double> > roots;

	if(fabs(a3) <= EPS) return Factor<2>(cs, EPS);
	a2 /= a3;
	a1 /= a3;
	a0 /= a3;

	double q = -(3 * a1 - a2 * a2) / 9;
	double r = -(9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
	double r2 = r * r;
	double q3 = q * q * q;

	if(r2 < q3) {
		double sqrQ = std::sqrt(q);
		double theta = std::acos(r / (sqrQ * q));
		double cTheta = std::cos(theta / 3) * sqrQ;
		double sTheta = std::sin(theta / 3) * sqrQ * SQRT_3 / 2;
		roots.push_back(-2 * cTheta);
		roots.push_back(-2 * (cTheta * 0.5 - sTheta));
		roots.push_back(-2 * (cTheta * 0.5 + sTheta));
	} else {
		double sqr = std::sqrt(r2 - q3);
		double t = -r + sqr;
		double s1 = t < 0 ? -pow(-t, 1.0 / 3) : pow(t, 1.0 / 3);
		t = -r - sqr;
		double s2 = t < 0 ? -pow(-t, 1.0 / 3) : pow(t, 1.0 / 3);
		roots.push_back(std::complex<double>(0, s1 + s2));
		s1 /= 2;
		s2 /= 2;
		roots.push_back(std::complex<double>(-s1 - s2, SQRT_3 * (s1 - s2)));
		roots.push_back(std::complex<double>(-s1 - s2, -SQRT_3 * (s1 - s2)));
	}
	roots[0].real(roots[0].real() - a2 / 3);
	roots[1].real(roots[1].real() - a2 / 3);
	roots[2].real(roots[2].real() - a2 / 3);

	return roots;
}
