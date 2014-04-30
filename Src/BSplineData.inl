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

/////////////////
// BSplineData //
/////////////////
// Support[i]:
//		Odd:  i +/- 0.5 * ( 1 + Degree )
//			i - 0.5 * ( 1 + Degree ) < 0
// <=>		i < 0.5 * ( 1 + Degree )
//			i + 0.5 * ( 1 + Degree ) > 0
// <=>		i > - 0.5 * ( 1 + Degree )
//			i + 0.5 * ( 1 + Degree ) > r
// <=>      i > r - 0.5 * ( 1 + Degree )
//			i - 0.5 * ( 1 + Degree ) < r
// <=>      i < r + 0.5 * ( 1 + Degree )
//		Even: i + 0.5 +/- 0.5 * ( 1 + Degree )
//			i - 0.5 * Degree < 0
// <=>		i < 0.5 * Degree
//			i + 1 + 0.5 * Degree > 0
// <=>		i > -1 - 0.5 * Degree
//			i + 1 + 0.5 * Degree > r
// <=>		i > r - 1 - 0.5 * Degree
//			i - 0.5 * Degree < r
// <=>		i < r + 0.5 * Degree

template<int Degree, class Real>
BSplineData<Degree, Real>::BSplineData():
	function_count_(0),
	sample_count_(0),
	value_tables_(NullPointer<Real>()) {
	SetBSplineElementIntegrals<Degree, Degree>(_vvIntegrals);
	SetBSplineElementIntegrals<Degree, Degree - 1>(_vdIntegrals);
	SetBSplineElementIntegrals<Degree - 1, Degree>(_dvIntegrals);
	SetBSplineElementIntegrals<Degree - 1, Degree - 1>(_ddIntegrals);
}

template<int Degree, class Real>
BSplineData<Degree, Real>::~BSplineData() {
	if(function_count_ && value_tables_) DeletePointer(value_tables_);
}

template<int Degree, class Real>
double BSplineData<Degree, Real>::Integrator::dot(int depth, int off1, int off2,
		bool d1, bool d2, bool childParent) const {
	if(depth < 0 || depth >= (int)iTables.size()) return 0.;
	typename Integrator::IntegralTables const& iTable = iTables[depth];
	int c;
	if(childParent) {
		c = off1 & 1;
		off1 >>= 1;
		--depth;
	}
	int d = off2 - off1;
	int res = 1 << depth;
	if(depth < 0 || off1 < 0 || off2 < 0 || off1 >= res || off2 >= res ||
			d < -Degree || d > Degree) return 0;
	int ii = off1 < Degree ? off1 : 
				off1 >= res - Degree ? 2 * Degree + off1 - (res - 1) :
				Degree;
	if(childParent)
		return (d1 && d2) ? iTable.dd_cpIntegrals[2 * ii + c][d + Degree] :
			d1 ? iTable.dv_cpIntegrals[2 * ii + c][d + Degree] :
			d2 ? iTable.vd_cpIntegrals[2 * ii + c][d + Degree] :
			iTable.vv_cpIntegrals[2 * ii + c][d + Degree];
	else
		return (d1 && d2) ? iTable.dd_ccIntegrals[ii][d + Degree] :
			d1 ? iTable.dv_ccIntegrals[ii][d + Degree] :
			d2 ? iTable.vd_ccIntegrals[ii][d + Degree] :
			iTable.vv_ccIntegrals[ii][d + Degree];
}

template<int Degree, class Real>
template<int Radius>
double BSplineData<Degree, Real>::CenterEvaluator<Radius>::value(int depth, int off1,
		int off2, bool d, bool childParent) const {
	if(depth < 0 || depth >= (int)vTables.size()) return 0.;
	int c = 1;
	if(childParent) {
		c = (off1 & 1) * 2;
		off1 >>= 1;
		--depth;
	}
	typename CenterEvaluator::ValueTables const& vTable = vTables[depth];
	int dd = off1 - off2;
	int res = 1 << depth;
	if(depth < 0 || off1 < 0 || off2 < 0 || off1 >= res || off2 >= res ||
			dd < -Radius || dd > Radius) return 0;
	int ii = off2 < Degree ? off2 :
		off2 >= res - Degree ? 2 * Degree + off2 - (res - 1) :
		Degree;
	return d ? vTable.dValues[ii][(dd + Radius) * 3 + c] :
			vTable.vValues[ii][(dd + Radius) * 3 + c];
}

template<int Degree, class Real>
template<int Radius>
double BSplineData<Degree, Real>::CornerEvaluator<Radius>::value(int depth, int off1,
		int c1, int off2, bool d, bool childParent) const {
	if(c1 < 0 || c1 >= 2) {
		std::cerr << "[WARNING] Clamping corner to {0, 1}" << std::endl;
		c1 = std::max(0, std::min(c1, 1));
	}
	if(depth < 0 || depth >= (int)vTables.size()) return 0.;
	int c = c1;
	if(childParent) {
		c = off1 & 1;
		off1 >>= 1;
		--depth;
	}
	typename CornerEvaluator::ValueTables const& vTable = vTables[depth];
	int dd = off1 - off2;
	int res = 1 << depth;
	if(depth < 0 || off1 < 0 || off2 < 0 || off1 >= res || off2 >= res ||
			dd < -Radius || dd > Radius) return 0;
	int ii = off2 < Degree ? off2 :
		off2 >= res - Degree ? 2 * Degree + off2 - res + 1 : Degree;
	return d ? vTable.dValues[ii][(dd + Radius) * 2 + c + c1] :
		vTable.vValues[ii][(dd + Radius) * 2 + c + c1];
}

template<int Degree, class Real>
void BSplineData<Degree, Real>::set(int maxDepth, BoundaryType boundaryType) {
	_boundaryType = boundaryType;

	depth_ = maxDepth;
	// [Warning] This assumes that the functions spacing is dual
	function_count_ = BinaryNode<double>::CumulativeCenterCount(depth_);
	sample_count_ = BinaryNode<double>::CenterCount(depth_) +
		BinaryNode<double>::CornerCount(depth_);
	baseFunctions = NewPointer<PPolynomial<Degree>>(function_count_);
	base_bsplines_ = NewPointer<BSplineComponents>(function_count_);

	BSplineComponents baseBSpline;
	PPolynomial<Degree> baseFunction = PPolynomial<Degree>::BSpline();
	for(int i = 0; i <= Degree; ++i)
		baseBSpline[i] = Polynomial<Degree>::BSplineComponent(i).shift(
				-(Degree + 1) / 2 + i - 0.5);
	StartingPolynomial<Degree> sPolys[Degree + 4];

	PPolynomial<Degree> leftBaseFunction;
	for(int i = 0; i < Degree + 3; ++i) {
		sPolys[i].start = -(Degree + 1) / 2 + i - 1.5;
		sPolys[i].p *= 0;
		if(i <= Degree) sPolys[i].p += baseBSpline[i].shift(-1) * (int)_boundaryType;
		if(i >= 1 && i <= Degree + 1) sPolys[i].p += baseBSpline[i - 1];
		for(int j = 0; j < i; ++j) sPolys[i].p -= sPolys[j].p;
	}
	leftBaseFunction.set(sPolys, Degree + 3);

	PPolynomial<Degree> rightBaseFunction;
	for(int i = 0; i < Degree + 3; ++i) {
		sPolys[i].start = -(Degree + 1) / 2 + i - 0.5;
		sPolys[i].p *= 0;
		if(i <= Degree) sPolys[i].p += baseBSpline[i];
		if(i >= 1 && i <= Degree + 1)
			sPolys[i].p += baseBSpline[i - 1].shift(1) * (int)_boundaryType;
		for(int j = 0; j < i; ++j) sPolys[i].p -= sPolys[j].p;
	}
	rightBaseFunction.set(sPolys, Degree + 3);

	PPolynomial<Degree> leftRightBaseFunction;
	for(int i = 0; i < Degree + 4; ++i) {
		sPolys[i].start = -(Degree + 1) / 2 + i - 1.5;
		sPolys[i].p *= 0;
		// The left-shifted B-spline
		if(i <= Degree) sPolys[i].p += baseBSpline[i].shift(-1) * (int)_boundaryType; 
		// The centered B-Spline
		if(i >= 1 && i <= Degree + 1) sPolys[i].p += baseBSpline[i - 1];
		// The right-shifted B-spline
		if(i >= 2 && i <= Degree + 2) sPolys[i].p += baseBSpline[i - 2].shift(1) *
			(int)_boundaryType; 
		for(int j = 0; j < i; ++j) sPolys[i].p -= sPolys[j].p;
	}
	leftRightBaseFunction.set(sPolys, Degree + 4);

	BSplineComponents leftBSpline;
	BSplineComponents rightBSpline;
	BSplineComponents leftRightBSpline;
	leftRightBSpline = leftBSpline = rightBSpline = baseBSpline;
	leftBSpline[1] += leftBSpline[2].shift(-1);
	leftBSpline[0] *= 0;
	rightBSpline[1] += rightBSpline[0].shift(1);
	rightBSpline[2] *= 0;
	leftRightBSpline[1] += leftRightBSpline[2].shift(-1) + leftRightBSpline[0].shift(1);
	leftRightBSpline[0] *= 0;
	leftRightBSpline[2] *= 0;

	for(size_t i = 0; i != function_count_; ++i) {
		auto caw = BinaryNode<double>::CenterAndWidth(i);
		baseFunctions[i] = baseFunction.scale(caw.second).shift(caw.first);
		base_bsplines_[i] = baseBSpline.scale(caw.second).shift(caw.first);
		if(_boundaryType != BoundaryType::None) {
			auto dao = BinaryNode<double>::DepthAndOffset(i);
			int r = 1 << dao.first;
			if(dao.second == 0 && dao.second == r - 1) {
				baseFunctions[i] = leftRightBaseFunction.scale(caw.second).shift(caw.first);
				base_bsplines_[i] = leftRightBSpline.scale(caw.second).shift(caw.first);
			} else if(dao.second == 0) {
				baseFunctions[i] = leftBaseFunction.scale(caw.second).shift(caw.first);
				base_bsplines_[i] = leftBSpline.scale(caw.second).shift(caw.first);
			} else if(dao.second == r - 1) {
				baseFunctions[i] = rightBaseFunction.scale(caw.second).shift(caw.first);
				base_bsplines_[i] = rightBSpline.scale(caw.second).shift(caw.first);
			}
		}
	}
}

template<int Degree, class Real>
template<bool D1, bool D2>
double BSplineData<Degree, Real>::_dot(int depth1, int off1, int depth2, int off2,
		bool inset) const {
	int const _Degree1 = D1 ? Degree - 1 : Degree;
	int const _Degree2 = D2 ? Degree - 1 : Degree;
	int sums[_Degree1 + 1][_Degree2 + 1]{ };

	int depth = std::max(depth1, depth2);

	BSplineElements<Degree> b1(1 << depth1, off1, _boundaryType,
			inset ? 1 << (depth1 - 2) : 0);
	BSplineElements<Degree> b2(1 << depth2, off2, _boundaryType,
			inset ? 1 << (depth2 - 2) : 0);

	BSplineElements<Degree> b;
	for(;depth1 < depth; ++depth1) {
		b = b1;
		b.upSample(b1);
	}
	for(;depth2 < depth; ++depth2) {
		b = b2;
		b.upSample(b2);
	}

	BSplineElements<Degree - 1> db1;
	BSplineElements<Degree - 1> db2;
	b1.differentiate(db1);
	b2.differentiate(db2);

	int start1 = -1;
	int end1 = -1;
	int start2 = -1;
	int end2 = -1;

	for(size_t i = 0; i != b1.size(); ++i)
		for(int j = 0; j <= Degree; ++j) {
			if(b1[i][j] && start1 == -1) start1 = i;
			if(b1[i][j]) end1 = i + 1;
			if(b2[i][j] && start2 == -1) start2 = i;
			if(b2[i][j]) end2 = i + 1;
		}

	if(start1 == end1 || start2 == end2 || start1 >= end2 || start2 >= end1) return 0;

	int start = std::max(start1, start2);
	int end = std::min(end1, end2);
	for(int i = start; i != end; ++i)
		for(int j = 0; j <= _Degree1; ++j)
			for(int k = 0; k <= _Degree2; ++k)
				sums[j][k] += (D1 ? db1[i][j] : b1[i][j]) * (D2 ? db2[i][k] : b2[i][k]);
	double _dot = 0;
	for(int j = 0; j <= _Degree1; ++j)
		for(int k = 0; k <= _Degree2; ++k)
			if(D1 && D2) _dot += _ddIntegrals[j][k] * sums[j][k];
			else if(D1) _dot += _dvIntegrals[j][k] * sums[j][k];
			else if(D2) _dot += _vdIntegrals[j][k] * sums[j][k];
			else _dot += _vvIntegrals[j][k] * sums[j][k];
	_dot /= b1.denominator() * b2.denominator();
	return D1 && D2 ? _dot * (1 << depth) :
		D1 || D2 ? _dot : _dot / (1 << depth);
}

template<int Degree, class Real>
double BSplineData<Degree, Real>::value(int depth, int off, double smoothingRadius,
		double s, bool d) const {
	if(off < 0 || off >= (1 << depth)) return 0;
	size_t idx = BinaryNode<Real>::CenterIndex(depth, off);

	PPolynomial<Degree + 1> function;
	if(smoothingRadius > 0) function = baseFunctions[idx].MovingAverage(smoothingRadius);
	else function = baseFunctions[idx];
	return d ? function.derivative()(s) : function(s);
}

template<int Degree, class Real>
void BSplineData<Degree, Real>::setIntegrator(Integrator& integrator, bool inset) const {
	integrator.iTables.resize(depth_ + 1);
	for(int d = 0; d <= depth_; ++d)
		for(int i = 0; i <= 2 * Degree; ++i)
			for(int j = -Degree; j <= Degree; ++j) {
				int res = 1 << d;
				int ii = i <= Degree ? i : i + res - 1 - 2 * Degree;
				integrator.iTables[d].vv_ccIntegrals[i][j + Degree] =
					_dot<false, false>(d, ii, d, ii + j, inset);
				integrator.iTables[d].dv_ccIntegrals[i][j + Degree] =
					_dot<true, false>(d, ii, d, ii + j, inset);
				integrator.iTables[d].vd_ccIntegrals[i][j + Degree] =
					_dot<false, true>(d, ii, d, ii + j, inset);
				integrator.iTables[d].dd_ccIntegrals[i][j + Degree] =
					_dot<true, true>(d, ii, d, ii + j, inset);
			}
	for(int d = 1; d <= depth_; ++d)
		for(int i = 0; i <= 2 * Degree; ++i)
			for(int j = -Degree; j <= Degree; ++j) {
				int res = 1 << d;
				int ii = i <= Degree ? i : i + res / 2 - 1 - 2 * Degree;
				for(int c = 0; c < 2; ++c) {
					integrator.iTables[d].vv_cpIntegrals[2 * i + c][j + Degree] =
						_dot<false, false>(d, 2 * ii + c, d - 1, ii + j, inset);
					integrator.iTables[d].dv_cpIntegrals[2 * i + c][j + Degree] =
						_dot<true, false>(d, 2 * ii + c, d - 1, ii + j, inset);
					integrator.iTables[d].vd_cpIntegrals[2 * i + c][j + Degree] =
						_dot<false, true>(d, 2 * ii + c, d - 1, ii + j, inset);
					integrator.iTables[d].dd_cpIntegrals[2 * i + c][j + Degree] =
						_dot<true, true>(d, 2 * ii + c, d - 1, ii + j, inset);
				}
			}
}

template<int Degree, class Real>
template<int Radius>
void BSplineData<Degree, Real>::setCenterEvaluator(CenterEvaluator<Radius>& evaluator,
		double smoothingRadius, double dSmoothingRadius) const {
	evaluator.vTables.resize(depth_ + 1);
	for(int d = 0; d <= depth_; ++d)
		for(int i = 0; i <= 2 * Degree; ++i)
			for(int j = -Radius; j <= Radius; ++j)
				for(int k = -1; k <= 1; ++k) {
					int res = 1 << d;
					int ii = i <= Degree ? i : i + res - 1 - 2 * Degree;
					double s = 0.5 + ii + j + 0.25 * k;
					evaluator.vTables[d].vValues[i][(j + Radius) * 3 + (k + 1)] =
						value(d, ii, smoothingRadius, s / res, false);
					evaluator.vTables[d].dValues[i][(j + Radius) * 3 + (k + 1)] =
						value(d, ii, dSmoothingRadius, s / res, true);
				}
}

template<int Degree, class Real>
template<int Radius>
void BSplineData<Degree, Real>::setCornerEvaluator(CornerEvaluator<Radius>& evaluator,
		double smoothingRadius, double dSmoothingRadius) const {
	evaluator.vTables.resize(depth_ + 1);
	for(int d = 0; d <= depth_; ++d)
		for(int i = 0; i <= 2 * Degree; ++i)
			for(int j = -Radius; j <= Radius; ++j)
				for(int k = 0; k <= 2; ++k) {
					int res = 1 << d;
					int ii = i <= Degree ? i : i + res - 1 - 2 * Degree;
					double s = ii + j + 0.5 * k;
					evaluator.vTables[d].vValues[i][(j + Radius) * 2 + k] =
						value(d, ii, smoothingRadius, s / res, false);
					evaluator.vTables[d].dValues[i][(j + Radius) * 2 + k] =
						value(d, ii, dSmoothingRadius, s / res, true);
				}
}

template< int Degree , class Real >
void BSplineData< Degree , Real >::setSampleSpan(int idx, int& start, int& end) const {
	auto dao = BinaryNode<double>::DepthAndOffset(idx);
	int res = 1 << dao.first;
	double _start = (dao.second + 0.5 - 0.5 * (Degree + 1)) / res;
	double _end   = (dao.second + 0.5 + 0.5 * (Degree + 1)) / res;
	// start / (sample_count_ - 1) > _start && (start - 1) / (sample_count_ - 1) <= _start
	// start > _start * (sample_count_ - 1) && start <= _start * (sample_count_ - 1) + 1
	// _start * (sample_count_ - 1) + 1 >= start > _start * (sample_count_ - 1)
	// end / (sample_count_ - 1) < _end && (end + 1) / (sample_count_ - 1) >= _end
	// end < _end * (sample_count_ - 1) && end >= _end * (sample_count_ - 1) - 1
	// _end * (sample_count_ - 1) > end >= _end * (sample_count_ - 1) - 1
	start = std::max<int>(floor(_start * (sample_count_ - 1) + 1), 0);
	end = std::min<int>(ceil(_end * (sample_count_ - 1) - 1), sample_count_ - 1);
}

template<int Degree, class Real>
void BSplineData<Degree, Real>::setValueTables() {
	if(value_tables_) DeletePointer(value_tables_);
	value_tables_ = NewPointer<Real>(function_count_ * sample_count_);
	for(size_t i = 0; i != function_count_; ++i)
		for(size_t j = 0; j != sample_count_; ++j)
			value_tables_[j * function_count_ + i] =
				baseFunctions[i]((double)j / (sample_count_ - 1));
}

/////////////////////
// BSplineElements //
/////////////////////
template<int Degree>
BSplineElements<Degree>::BSplineElements(int res, int offset, BoundaryType boundary,
		int inset):
	denominator_(1),
	elems(res) {

	for(int i = 0; i <= Degree; ++i) {
		int idx = -_off + offset + i;
		if(idx >= 0 && idx < res) elems[idx][i] = 1;
	}
	if(boundary != BoundaryType::None) {
		_addLeft(offset - 2 * res, boundary);
		_addRight(offset + 2 * res, boundary);
		if(Degree & 1) {
			_addLeft(offset - res, boundary);
			_addRight(offset + res, boundary);
		} else {
			_addLeft(-offset - 1, boundary);
			_addRight(-offset - 1 + 2 * res, boundary);
		}
	}
	if(inset)
		for(int i = 0; i < inset && i < res; ++i)
			for(int j = 0; j <= Degree; ++j)
				elems[i][j] = elems[res - 1 - i][j] = 0;
}

template<int Degree>
void BSplineElements<Degree>::_addLeft(int offset, BoundaryType boundary ) {
	int res = elems.size();
	bool set = false;
	for(int i = 0; i <= Degree; ++i) {
		int idx = -_off + offset + i;
		if(idx >= 0 && idx < res) {
			elems[idx][i] += (int)boundary;
			set = true;
		}
	}
	if(set) _addLeft(offset - 2 * res, boundary);
}

template<int Degree>
void BSplineElements<Degree>::_addRight(int offset, BoundaryType boundary) {
	int res = elems.size();
	bool set = false;
	for(int i = 0; i <= Degree; ++i) {
		int idx = -_off + offset + i;
		if(idx >= 0 && idx < res) {
			elems[idx][i] += (int)boundary;
			set = true;
		}
	}
	if(set) _addRight(offset + 2 * res, boundary);
}

template<int Degree>
void BSplineElements<Degree>::upSample(BSplineElements<Degree>& high) const {
	std::cerr << "[ERROR] B-spline up-sampling not supported for degree "
		<< Degree << std::endl;
	exit(0);
}

template<>
void BSplineElements<1>::upSample(BSplineElements<1>& high) const {
	high.elems.resize(elems.size() * 2);
	high.elems.assign(high.elems.size(), BSplineElementCoefficients<1>());
	for(size_t i = 0; i != elems.size(); ++i) {
		high.elems[2 * i][0] += elems[i][0] + elems[i][1];
		high.elems[2 * i][1] += 2 * elems[i][1];
		high.elems[2 * i + 1][0] += 2 * elems[i][0];
		high.elems[2 * i + 1][1] += elems[i][0] + elems[i][1];
	}
	high.denominator_ = denominator_ * 2;
}

template<>
void BSplineElements<2>::upSample(BSplineElements<2>& high) const {
	/*
	    /----\
	   /      \
	  /        \  = 1  /--\       +3    /--\     +3      /--\   +1        /--\
	 /          \     /    \           /    \           /    \           /    \
	 |----------|     |----------|   |----------|   |----------|   |----------|
	*/
	high.elems.resize(elems.size() * 2);
	high.elems.assign(high.elems.size(), BSplineElementCoefficients<2>());
	for(size_t i = 0; i != elems.size(); ++i) {
		high.elems[2 * i][0] += elems[i][0] + 3 * elems[i][1];
		high.elems[2 * i][1] += 3 * elems[i][1] + elems[i][2];
		high.elems[2 * i][2] += elems[i][1] + 3 * elems[i][2];
		high.elems[2 * i + 1][0] += 3 * elems[i][0] + elems[i][1];
		high.elems[2 * i + 1][1] += 1 * elems[i][0] + 3 * elems[i][1];
		high.elems[2 * i + 1][2] += 3 * elems[i][1] + 1 * elems[i][2];
	}
	high.denominator_ = denominator_ * 4;
}

template<int Degree>
void BSplineElements<Degree>::differentiate(BSplineElements<Degree - 1>& d) const {
	d.elems.resize(elems.size());
	d.elems.assign(d.elems.size(), BSplineElementCoefficients<Degree - 1>());
	for(size_t i = 0; i != elems.size(); ++i)
		for(int j = 0; j <= Degree; ++j) {
			if(j - 1 >= 0) d.elems[i][j - 1] -= elems[i][j];
			if(j < Degree) d.elems[i][j] += elems[i][j];
	}
	d.denominator_ = denominator_;
}

template<int Degree>
void BSplineElements<Degree>::print() const {
	for(int i = 0; i != elems.size(); ++i) {
		std::cout << i << "]";
		for(int j = 0; j <= Degree; ++j) std::cout << " " << elems[i][j];
		std::cout << " (" << denominator_ << ")" << std::endl;
	}
}

// If we were really good, we would implement this integral table to store
// rational values to improve precision...
template<int Degree1, int Degree2>
void SetBSplineElementIntegrals(double integrals[Degree1 + 1][Degree2 + 1]) {
	for(int i = 0; i <= Degree1; ++i) {
		Polynomial<Degree1> p1 = Polynomial<Degree1>::BSplineComponent(i);
		for(int j = 0; j <= Degree2; ++j) {
			Polynomial<Degree2> p2 = Polynomial<Degree2>::BSplineComponent(j);
			integrals[i][j] = (p1 * p2).integral(0, 1);
		}
	}
}
