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

#include "Vector.h"
#include "Array.h"

//#define NEW_MATRIX_CODE

template<class T>
struct MatrixEntry {
	MatrixEntry(): N(-1), Value(0) { }
	MatrixEntry(int i, T v = T()): N(i), Value(v) { }
	int N;
	T Value;
};

template<class T>
class MapReduceVector {
public:
	MapReduceVector(): _dim(0) { }
	~MapReduceVector() {
		if(_dim) for(size_t t = 0; t != out.size(); ++t) delete[] out[t];
		out.resize(0);
	}
	T* operator[](int t) { return out[t]; }
	T const* operator[](int t) const { return out[t]; }
	int threads() const { return out.size(); }
	void resize(size_t threads, int dim) {
		if(threads != out.size() || _dim < dim) {
			for(size_t t = 0; t != out.size(); ++t) delete[] out[t];
			out.resize(threads);
			for(size_t t = 0; t != out.size(); ++t) out[t] = new T[dim];
			_dim = dim;
		}
	}
private:
	std::vector<T*> out;
	int _dim;
};

template<class T>
class SparseSymmetricMatrix {
public:
	int Rows() const { return rows; }
	int Columns() const { return _maxEntriesPerRow; }

	Pointer(MatrixEntry<T>) operator[](int idx) { return m_ppElements[idx]; }
	ConstPointer(MatrixEntry<T>) operator[](int idx) const { return m_ppElements[idx]; }

	int& rowSize(int i) { return rowSizes_[i]; }

	SparseSymmetricMatrix();
	SparseSymmetricMatrix(SparseSymmetricMatrix const& M);

	~SparseSymmetricMatrix();

	SparseSymmetricMatrix<T>& operator=(SparseSymmetricMatrix<T> M) { swap(M); return *this; }

	void swap(SparseSymmetricMatrix<T>& M);

	void Resize(int rows);
	void SetRowSize(int row, int count);
	int Entries() const;

	template<class T2>
	Vector<T2> operator*(Vector<T2> const& V) const;

#ifdef NEW_MATRIX_CODE
	template<class T2>
	static int Solve(SparseSymmetricMatrix<T> const& M, Vector<T2> const& b, int iters, Vector<T2>& solution,
			T2 eps, bool reset, int threads, bool addDCTerm);
#else
	template<class T2>
	static int Solve(SparseSymmetricMatrix<T> const& M, Vector<T2> const& b, int iters, Vector<T2>& solution,
			MapReduceVector<T2>& scratch, T2 eps, bool reset, bool addDCTerm);
#endif
private:
#ifdef NEW_MATRIX_CODE
	template<class T2>
	void Multiply(Vector<T2> const& In, Vector<T2>& Out, bool addDCTerm, int threads) const;
#else
	template<class T2>
	void Multiply(Vector<T2> const& In, Vector<T2>& Out, MapReduceVector<T2>& OutScratch, bool addDCTerm) const;
#endif
private:
	Pointer(int) rowSizes_;
	Pointer(Pointer(MatrixEntry<T>)) m_ppElements;
	bool _contiguous;
	int _maxEntriesPerRow;
	int rows;
};

#include "SparseMatrix.inl"
