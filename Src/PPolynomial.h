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

#include <vector>

#include "Polynomial.h"

template<int Degree>
struct StartingPolynomial {
	Polynomial<Degree> p;
	double start;

	StartingPolynomial(): start(0) { }

	StartingPolynomial(Polynomial<Degree> const& p, double start):
		p(p), start(start) { }

	template<int Degree2>
	StartingPolynomial(StartingPolynomial<Degree2> const& p):
		p(p.p), start(p.start) { }

	StartingPolynomial scale(double s) const
		{ return StartingPolynomial(p.scale(s), start * s); }
	StartingPolynomial shift(double t) const
		{ return StartingPolynomial(p.shift(t), start + t); }
};

template<int D>
bool operator<(StartingPolynomial<D> const& p1, StartingPolynomial<D> const& p2)
	{ return p1.start < p2.start; }

template<int Degree>
class PPolynomial {
public:
	static PPolynomial BSpline(double radius = 0.5);
public:
	PPolynomial() { }
	// Note: this constructor will sort the elements in sps
	PPolynomial(StartingPolynomial<Degree>* sps, int count);

	template<int Degree2>
	PPolynomial(PPolynomial<Degree2> const&);

	template<int Degree2>
	PPolynomial& operator=(PPolynomial<Degree2> const& p);

	void swap(PPolynomial& p);

	double operator()(double) const;

	PPolynomial scale(double) const;
	PPolynomial shift(double) const;

	PPolynomial<Degree - 1> derivative() const;

	PPolynomial<Degree + 1> MovingAverage(double radius) const;
private:
	explicit PPolynomial(size_t size);

	PPolynomial& operator/=(double);

	friend PPolynomial operator/(PPolynomial p, double s) { return p /= s; }
private:
	template<int Degree2>
	friend class PPolynomial;

	std::vector<StartingPolynomial<Degree> > polys_;
};

#include "PPolynomial.inl"
