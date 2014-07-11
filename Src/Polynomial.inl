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

#include <float.h>
#include <math.h>
#include <algorithm>

#include "Factor.h"

////////////////
// Polynomial //
////////////////

template<int Degree>
Polynomial<Degree>::Polynomial(double cs[Degree + 1]) {
	std::copy(cs, cs + Degree + 1, coefficients);
}

template<int Degree>
template<int Degree2>
Polynomial<Degree>::Polynomial(Polynomial<Degree2> const& P) {
	int i;
	for(i = 0; i <= Degree && i <= Degree2; ++i)
		coefficients[i] = P[i];
	for(; i <= Degree; ++i)
		coefficients[i] = 0;
}

template<int Degree>
void Polynomial<Degree>::swap(Polynomial& p) {
	using std::swap;
	swap(coefficients, p.coefficients);
}

template<int Degree>
template<int Degree2>
Polynomial<Degree>& Polynomial<Degree>::operator=(Polynomial<Degree2> const& p) {
	Polynomial<Degree> tmp(p);
	swap(tmp);
	return *this;
}

template<int Degree>
Polynomial<Degree - 1> Polynomial<Degree>::derivative() const {
	Polynomial<Degree - 1> p;
	for(int i = 0; i != Degree; ++i) p[i] = coefficients[i + 1] * (i + 1);
	return p;
}

template<int Degree>
Polynomial<Degree + 1> Polynomial<Degree>::integral() const {
	Polynomial<Degree + 1> p;
	p[0] = 0;
	for(int i = 0; i <= Degree; ++i) p[i + 1] = coefficients[i] / (i + 1);
	return p;
}

template<int Degree>
double Polynomial<Degree>::operator()(double t) const {
	double v = 0;
	for(int d = Degree; d >= 0; --d) v = v * t + coefficients[d];
	return v;
}

// TODO: Maybe just integral()(tMax) - integral()(tMin)
template<int Degree>
double Polynomial<Degree>::integral(double tMin, double tMax) const {
	double t1 = tMin;
	double t2 = tMax;
	double v = 0;
	for(int i = 0; i <= Degree; ++i) {
		v += coefficients[i] * (t2 - t1) / (i + 1);
		if(t1 != -DBL_MAX && t1 != DBL_MAX) t1 *= tMin;
		if(t2 != -DBL_MAX && t2 != DBL_MAX) t2 *= tMax;
	}
	return v;
}

template<int Degree>
Polynomial<Degree>& Polynomial<Degree>::operator+=(Polynomial<Degree> const& p) {
	for(int i = 0; i <= Degree; ++i) coefficients[i] += p.coefficients[i];
	return *this;
}

template<int Degree>
Polynomial<Degree>& Polynomial<Degree>::operator-=(Polynomial<Degree> const& p) {
	for(int i = 0; i <= Degree; ++i) coefficients[i] -= p.coefficients[i];
	return *this;
}

template<int Degree>
Polynomial<Degree> Polynomial<Degree>::operator-() const {
	Polynomial q;
	for(int i = 0; i <= Degree; ++i) q.coefficients[i] = -coefficients[i];
	return q;
}

template<int Degree>
template<int Degree2>
Polynomial<Degree + Degree2> Polynomial<Degree>::operator*(Polynomial<Degree2> const& p) const {
	Polynomial<Degree + Degree2> q;
	for(int i = 0; i <= Degree; ++i)
		for(int j = 0; j <= Degree2; ++j)
			q[i + j] += coefficients[i] * p[j];
	return q;
}

template<int Degree>
Polynomial<Degree>& Polynomial<Degree>::operator+=(double s) {
	coefficients[0] += s;
	return *this;
}

template<int Degree>
Polynomial<Degree>& Polynomial<Degree>::operator-=(double s) {
	coefficients[0] -= s;
	return *this;
}

template<int Degree>
Polynomial<Degree>& Polynomial<Degree>::operator*=(double s) {
	for(int i = 0; i <= Degree; ++i) coefficients[i] *= s;
	return *this;
}

template<int Degree>
Polynomial<Degree>& Polynomial<Degree>::operator/=(double s) {
	for(int i = 0; i <= Degree; ++i) coefficients[i] /= s;
	return *this;
}

template<int Degree>
Polynomial<Degree> Polynomial<Degree>::scale(double s) const {
	Polynomial q = *this;
	double s2 = 1.0;
	for(int i = 0; i <= Degree; ++i) {
		q.coefficients[i] *= s2;
		s2 /= s;
	}
	return q;
}

template<int Degree>
Polynomial<Degree> Polynomial<Degree>::shift(double t) const {
	Polynomial<Degree> q;
	for(int i = 0; i <= Degree; ++i) {
		double temp = 1;
		for(int j = i; j >= 0; j--) {
			q.coefficients[j] += coefficients[i] * temp;
			temp *= -t * j;
			temp /= i - j + 1;
		}
	}
	return q;
}

template<int Degree>
std::vector<double> Polynomial<Degree>::getSolutions(double c, double EPS) const {
	double cs[Degree + 1];
	std::copy(coefficients, coefficients + Degree + 1, cs);
	cs[0] -= c;
	std::vector<std::complex<double> > rs = Factor<Degree>(cs, EPS);
	std::vector<double> roots;
	for(std::vector<std::complex<double> >::iterator r = rs.begin(); r != rs.end(); ++r)
		if(std::abs(r->imag()) <= EPS) roots.push_back(r->real());
	return roots;
}

template<>
Polynomial<0> Polynomial<0>::BSplineComponent(int) {
	Polynomial p;
	p[0] = 1;
	return p;
}

template<int Degree>
Polynomial<Degree> Polynomial<Degree>::BSplineComponent(int i) {
	Polynomial p;
	if(i > 0) {
		Polynomial<Degree> _p = Polynomial<Degree - 1>::BSplineComponent(i - 1).integral();
		p -= _p;
		p[0] += _p(1);
	}
	if(i < Degree) {
		Polynomial<Degree> _p = Polynomial<Degree - 1>::BSplineComponent(i).integral();
		p += _p;
	}
	return p;
}
