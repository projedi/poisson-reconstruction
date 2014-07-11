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

#pragma once

#include <algorithm>
#include <vector>

template<int Degree>
class Polynomial {
public:
	Polynomial() { memset(coefficients, 0, sizeof(double) * (Degree + 1)); }

	Polynomial(double cs[Degree + 1]);

	template<int Degree2>
	Polynomial(Polynomial<Degree2> const& P);

	template<int Degree2>
	Polynomial& operator=(Polynomial<Degree2> const& p);

	void swap(Polynomial& p);

	double& operator[](size_t i) { return coefficients[i]; }
	double const& operator[](size_t i) const { return coefficients[i]; }

	double operator()(double t) const;
	double integral(double tMin, double tMax) const;

	bool operator==(Polynomial const& p) const
		{ return std::equal(coefficients.begin(), coefficients.end(), p.coefficients.begin()); }
	bool operator!=(Polynomial const& p) const { return !(*this == p); }

	bool isZero() const
		{ return std::all_of(coefficients, coefficients + Degree + 1, [](double p) { return p == 0; }); }
	void setZero() { std::fill(coefficients.begin(), coefficients.end(), 0); }

	Polynomial& operator+=(Polynomial const& p);
	Polynomial& operator-=(Polynomial const& p);
	Polynomial operator-() const;

	template<int Degree2>
	Polynomial<Degree + Degree2> operator*(Polynomial<Degree2> const& p2) const;

	Polynomial& operator+=(double s);
	Polynomial& operator-=(double s);
	Polynomial& operator*=(double s);
	Polynomial& operator/=(double s);

	Polynomial scale(double s) const;
	Polynomial shift(double t) const;

	Polynomial<Degree - 1> derivative() const;
	Polynomial<Degree + 1> integral() const;

	std::vector<double> getSolutions(double c, double EPS) const;

	static Polynomial BSplineComponent(int i);
private:
	double coefficients[Degree + 1];
};

template<int Degree>
Polynomial<Degree> operator+(Polynomial<Degree> p1, Polynomial<Degree> const& p2) { return p1 += p2; }

template<int Degree>
Polynomial<Degree> operator-(Polynomial<Degree> p1, Polynomial<Degree> const& p2) { return p1 -= p2; }

template<int Degree>
Polynomial<Degree> operator+(Polynomial<Degree> p, double s) { return p += s; }

template<int Degree>
Polynomial<Degree> operator-(Polynomial<Degree> p, double s) { return p -= s; }

template<int Degree>
Polynomial<Degree> operator*(Polynomial<Degree> p, double s) { return p *= s; }

template<int Degree>
Polynomial<Degree> operator/(Polynomial<Degree> p, double s) { return p /= s; }

#include "Polynomial.inl"
