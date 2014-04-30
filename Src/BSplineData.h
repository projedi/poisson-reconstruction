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

#include "Array.h"
#include "PPolynomial.h"

template<int Degree>
struct BSplineElementCoefficients {
	BSplineElementCoefficients(): coeffs{} { }
	int operator[](int idx) const { return coeffs[idx]; }
	int& operator[](int idx) { return coeffs[idx]; }
private:
	int coeffs[Degree + 1];
};

enum class BoundaryType {
	None = 0,
	Dirichlet = -1,
	Neumann = 1
};

template<int Degree>
struct BSplineElements {
	BSplineElements(): denominator_(1) { }
	BSplineElements(int res, int offset, BoundaryType boundary = BoundaryType::None,
			int inset = 0);

	void upSample(BSplineElements& high) const;
	void differentiate(BSplineElements<Degree - 1>& d) const;

	void print() const;

	int denominator() const { return denominator_; }

	size_t size() const { return elems.size(); }
	BSplineElementCoefficients<Degree>& operator[](size_t i) { return elems[i]; }
	BSplineElementCoefficients<Degree> const& operator[](size_t i) const
		{ return elems[i]; }
private:
	void _addLeft(int offset, BoundaryType boundary);
	void _addRight(int offset, BoundaryType boundary);
private:
	friend struct BSplineElements<Degree + 1>;
	static int const _off = (Degree + 1) / 2;
	// Coefficients are ordered as "/" "-" "\"
	int denominator_;
	std::vector<BSplineElementCoefficients<Degree>> elems;
};

template<int Degree, class Real>
class BSplineData {
public:
	BSplineData();
	~BSplineData();

	struct Integrator {
		double dot(int depth, int off1, int off2, bool d1, bool d2,
				bool childParent = false ) const;
	private:
		friend class BSplineData<Degree, Real>;
		struct IntegralTables {
			double vv_ccIntegrals[2 * Degree + 1][2 * Degree + 1];
			double dv_ccIntegrals[2 * Degree + 1][2 * Degree + 1];
			double vd_ccIntegrals[2 * Degree + 1][2 * Degree + 1];
			double dd_ccIntegrals[2 * Degree + 1][2 * Degree + 1];
			double vv_cpIntegrals[(2 * Degree + 1) * 2][2 * Degree + 1];
			double dv_cpIntegrals[(2 * Degree + 1) * 2][2 * Degree + 1];
			double vd_cpIntegrals[(2 * Degree + 1) * 2][2 * Degree + 1];
			double dd_cpIntegrals[(2 * Degree + 1) * 2][2 * Degree + 1];
		};
		std::vector<IntegralTables> iTables;
	};
	void setIntegrator(Integrator& integrator, bool inset) const;

	template<int Radius>
	struct CenterEvaluator {
		double value(int depth, int off1, int off2, bool d, bool childParent = false) const;
	private:
		friend class BSplineData<Degree, Real>;
		struct ValueTables {
			double vValues[2 * Degree + 1][3 * (2 * Radius + 1)];
			double dValues[2 * Degree + 1][3 * (2 * Radius + 1)];
		};
		std::vector<ValueTables> vTables;
	};
	template<int Radius>
	void setCenterEvaluator(CenterEvaluator<Radius>& evaluator, double smoothingRadius,
			double dSmoothingRadius) const;

	template<int Radius>
	struct CornerEvaluator {
		double value(int depth, int off1, int c1, int off2, bool d,
				bool childParent = false) const;
	private:
		friend class BSplineData<Degree, Real>;
		struct ValueTables {
			double vValues[2 * Degree + 1][4 * Radius + 3];
			double dValues[2 * Degree + 1][4 * Radius + 3];
		};
		std::vector<ValueTables> vTables;
	};
	template<int Radius>
	void setCornerEvaluator(CornerEvaluator<Radius>& evaluator, double smoothingRadius,
			double dSmoothingRadius) const;

	void setValueTables();

	void setSampleSpan(int idx, int& start, int& end) const;

	/********************************************************
	 * Sets the translates and scales of the basis function
	 * up to the prescribed depth
	 * <maxDepth> the maximum depth
	 ********************************************************/
	void set(int maxDepth, BoundaryType boundaryType);

	int depth() const { return depth_; }
	size_t functionCount() const { return function_count_; }
	Real& valueTables(size_t idx) { return value_tables_[idx]; }
	Real const& valueTables(size_t idx) const { return value_tables_[idx]; }
	Polynomial<Degree>& baseBSplines(size_t i, size_t j) { return base_bsplines_[i][j]; }
	Polynomial<Degree> const& baseBSplines(size_t i, size_t j) const
		{ return base_bsplines_[i][j]; }
private:
	template<bool D1, bool D2>
	double _dot(int depth1, int off1, int depth2, int off2, bool inset) const;

	double value(int depth, int off, double smoothingRadius, double s, bool d) const;
private:
	struct BSplineComponents {
		Polynomial<Degree> polys[Degree + 1];
		Polynomial<Degree>& operator[](int idx) { return polys[idx]; }
		Polynomial<Degree> const& operator[](int idx) const { return polys[idx]; }
		void printnl() const { for(int d = 0; d <= Degree; ++d) polys[d].printnl(); }
		BSplineComponents scale(double s) const {
			BSplineComponents b;
			for(int d = 0; d <= Degree; ++d)
				b[d] = polys[d].scale(s);
			return b;
		}
		BSplineComponents shift(double s) const {
			BSplineComponents b;
			for(int d = 0; d <= Degree; ++d)
				b[d] = polys[d].shift(s);
			return b;
		}
	};

	BoundaryType _boundaryType;
	double _vvIntegrals[Degree + 1][Degree + 1];
	double _vdIntegrals[Degree + 1][Degree];
	double _dvIntegrals[Degree][Degree + 1];
	double _ddIntegrals[Degree][Degree];

	int depth_;
	size_t function_count_;
	size_t sample_count_;
	Pointer(Real) value_tables_;
	Pointer(PPolynomial<Degree>) baseFunctions;
	Pointer(BSplineComponents) base_bsplines_;
};

template<int Degree1, int Degree2>
void SetBSplineElementIntegrals(double integrals[Degree1 + 1][Degree2 + 1]);

#include "BSplineData.inl"
