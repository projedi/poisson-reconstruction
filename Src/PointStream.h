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
prior writften permission. 

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

#include <fstream>

template<class Real>
class PointStream {
public:
	virtual ~PointStream() { }
	virtual void reset() = 0;
	virtual bool nextPoint(Point3D<Real>& p, Point3D<Real>& n) = 0;
};

template<class Real>
class ASCIIPointStream: public PointStream<Real> {
public:
	ASCIIPointStream(std::string const& filename);
	void reset() override { file_.clear(); file_.seekg(0); }
	bool nextPoint(Point3D<Real>& p, Point3D<Real>& n) override;
private:
	std::ifstream file_;
};

// TODO: Implementation is untouched because there is no test data
template<class Real>
class BinaryPointStream: public PointStream<Real> {
public:
	BinaryPointStream(std::string const& filename);
	~BinaryPointStream();
	void reset() override;
	bool nextPoint(Point3D<Real>& p, Point3D<Real>& n) override;
private:
	FILE* _fp;
	static int const POINT_BUFFER_SIZE = 1024;
	Real _pointBuffer[POINT_BUFFER_SIZE * 2 * 3];
	int _pointsInBuffer;
	int _currentPointIndex;
};

template<class Real>
class PLYPointStream: public PointStream<Real> {
public:
	PLYPointStream(std::string const& filename): _fileName(filename), _ply(nullptr) { reset(); }
	~PLYPointStream() { _free(); }
	void reset() override;
	bool nextPoint(Point3D<Real>& p, Point3D<Real>& n) override;
private:
	void _free( void );
private:
	std::string _fileName;
	PlyFile* _ply;
	int _nr_elems;
	char **_elist;

	int _pCount;
	int _pIdx;
};

#include "PointStream.inl"
