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

#include <numeric>

#include "Vector.h"

template<class T>
struct MatrixEntry {
	MatrixEntry(): N(-1), Value(0) { }
	MatrixEntry(int i, T v = T()): N(i), Value(v) { }
	int N;
	T Value;
};

template<class T>
class SparseSymmetricMatrix {
public:
	int Rows() const { return m_ppElements.size(); }

	MatrixEntry<T> const& at(int i, int j) const {
		if(j >= rowSizes_[i]) {
			//std::cerr << "[WARNING] accessing to the right of rowSize" << std::endl;
		}
		if((size_t)j >= m_ppElements[i].size()) {
			std::cerr << "[ERROR] accessing to the right of m_ppElements[i]" << std::endl;
			throw "OHAI";
		}
		return m_ppElements[i][j];
	}
	MatrixEntry<T>& at(int i, int j)
		{ return const_cast<MatrixEntry<T>&>(static_cast<SparseSymmetricMatrix const&>(*this).at(i, j)); }

	int& rowSize(int i) { return rowSizes_[i]; }

	void Resize(int rows) { m_ppElements.resize(rows); rowSizes_.resize(rows); }
	void SetRowSize(int row, int count) { m_ppElements[row].resize(count); }
	int Entries() const { return std::accumulate(rowSizes_.begin(), rowSizes_.end(), 0); }

	template<class T2>
	Vector<T2> operator*(Vector<T2> const& V) const;

	T Norm(size_t Ln) const;

	template<class T2>
	static int Solve(SparseSymmetricMatrix<T> const& M, Vector<T2> const& b, int iters, Vector<T2>& solution,
			T2 eps, bool reset, int threads, bool addDCTerm);
private:
	template<class T2>
	void Multiply(Vector<T2> const& In, Vector<T2>& Out, bool addDCTerm, int threads) const;
private:
	std::vector<int> rowSizes_;
	std::vector<std::vector<MatrixEntry<T> > > m_ppElements;
};

#include "SparseMatrix.inl"
