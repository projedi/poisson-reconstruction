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

template<int Degree>
PPolynomial<Degree>::PPolynomial(size_t size): polys_(size) { }

template<int Degree>
template<int Degree2>
PPolynomial<Degree>::PPolynomial(PPolynomial<Degree2> const& p):
	polys_(p.polys_.begin(), p.polys_.end()) { }

template<int Degree>
PPolynomial<Degree>::PPolynomial(StartingPolynomial<Degree>* sps, int count):
	PPolynomial(count) {
	std::sort(sps, sps + count);
	int c = 0;
	for(int i = 0; i != count; ++i) {
		if(!c || sps[i].start != polys_[c - 1].start) polys_[c++] = sps[i];
		else polys_[c - 1].p += sps[i].p;
	}
	polys_.resize(c);
}

template<int Degree>
template<int Degree2>
PPolynomial<Degree>& PPolynomial<Degree>::operator=(PPolynomial<Degree2> const& p) {
	PPolynomial tmp(p);
	swap(tmp);
	return *this;
}

template<int Degree>
double PPolynomial<Degree>::operator()(double t) const {
	double v = 0;
	for(size_t i = 0; i != polys_.size() && t > polys_[i].start; ++i) v += polys_[i].p(t);
	return v;
}

template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::scale(double s) const {
	PPolynomial q(polys_.size());
	for(size_t i = 0; i != polys_.size(); ++i) q.polys_[i] = polys_[i].scale(s);
	return q;
}

template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::shift(double s) const {
	PPolynomial q(polys_.size());
	for(size_t i = 0; i != polys_.size(); ++i) q.polys_[i] = polys_[i].shift(s);
	return q;
}

template<int Degree>
PPolynomial<Degree - 1> PPolynomial<Degree>::derivative() const {
	PPolynomial<Degree - 1> q(polys_.size());
	for(size_t i = 0; i != polys_.size(); ++i)
		q.polys_[i] = StartingPolynomial<Degree - 1>(polys_[i].p.derivative(),
				polys_[i].start);
	return q;
}

template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator/=(double s) {
	for(size_t i = 0; i != polys_.size(); ++i) polys_[i].p /= s;
	return *this;
}

template<>
PPolynomial<0> PPolynomial<0>::BSpline(double radius) {
	PPolynomial q(2);
	q.polys_[0].start = -radius;
	q.polys_[0].p.coefficients[0] = 1.0;
	q.polys_[1].start = radius;
	q.polys_[1].p.coefficients[0] = -1.0;
	return q;
}

template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::BSpline(double radius) {
	return PPolynomial<Degree - 1>::BSpline().MovingAverage(radius);
}

template<int Degree>
PPolynomial<Degree + 1> PPolynomial<Degree>::MovingAverage(double radius) const {

	std::vector<StartingPolynomial<Degree + 1>> sps(polys_.size() * 2);

	for(size_t i = 0; i != polys_.size(); ++i) {
		Polynomial<Degree + 1> p = polys_[i].p.integral() -
			polys_[i].p.integral()(polys_[i].start);
		sps[2 * i].start = polys_[i].start - radius;
		sps[2 * i + 1].start = polys_[i].start + radius;
		sps[2 * i].p = p.shift(-radius);
		sps[2 * i + 1].p = -p.shift(radius);
	}
	PPolynomial<Degree + 1> A(sps.data(), polys_.size() * 2);
	return A / (2 * radius);
}

template<int Degree>
void PPolynomial<Degree>::swap(PPolynomial& p) {
	using std::swap;
	swap(polys_, p.polys_);
}
