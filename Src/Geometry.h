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

#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>

template<class Real>
struct Point3D {
	Point3D(): coords{} { }
	Point3D(Real v0, Real v1, Real v2): coords{v0, v1, v2} { }
	template<class Real2>
	explicit Point3D(Point3D<Real2> const& p):
		coords{(Real)p[0], (Real)p[1], (Real)p[2]} { }

	Real& operator[](int i) { return coords[i]; }
	Real const& operator[](int i) const { return coords[i]; }

// TODO: in PlyFile offsetof is used on this field. MUST KILL WITH FIRE.
public:
	Real coords[3];
};

template<class Real>
Point3D<Real> operator-(Point3D<Real> const& p)
	{ return Point3D<Real>(-p[0], -p[1], -p[2]); }

template<class Real>
Point3D<Real>& operator+=(Point3D<Real>& p1, Point3D<Real> const& p2);

template<class Real>
Point3D<Real>& operator*=(Point3D<Real>& p, Real r);

template<class Real>
Point3D<Real>& operator/=(Point3D<Real>& p, Real r);

template<class Real>
Point3D<Real>& operator-=(Point3D<Real>& p1, Point3D<Real> const& p2)
	{ return p1 += -p2; }

template<class Real>
Point3D<Real> operator+(Point3D<Real> p1, Point3D<Real> const& p2) { return p1 += p2; }

template<class Real>
Point3D<Real> operator-(Point3D<Real> p1, Point3D<Real> const& p2) { return p1 -= p2; }

template<class Real1, class Real2>
Point3D<Real1> operator*(Point3D<Real1> p, Real2 r) { return p *= (Real1)r; }

template<class Real1, class Real2>
Point3D<Real1> operator/(Point3D<Real1> const& p, Real2 r) { return p * (1.0 / r); }

template<class Real>
double Dot(Point3D<Real> const& p1, Point3D<Real> const& p2)
	{ return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]; }

template<class Real>
double SquareLength(Point3D<Real> const& p) { return Dot(p, p); }

template<class Real>
double Length(Point3D<Real> const& p) { return std::sqrt(SquareLength(p)); }

template<class Real>
double SquareDistance(Point3D<Real> const& p1, Point3D<Real> const& p2)
	{ return SquareLength(p1 - p2); }

template<class Real>
double Distance(Point3D<Real> const& p1, Point3D<Real> const& p2)
	{ return std::sqrt(SquareDistance(p1, p2)); }

template <class Real>
Point3D<Real> CrossProduct(Point3D<Real> const& p1, Point3D<Real> const& p2);

template<class Real>
double TriangleArea(Point3D<Real> const& v1, Point3D<Real> const& v2, Point3D<Real> const& v3) {
	return Length(CrossProduct(v2 - v1, v3 - v1)) / 2;
}

template<class Real, int Dim>
struct XForm {
	XForm(): coords_{} { }

	static XForm Identity();

	Real& operator()(int i, int j) { return coords_[i][j]; }
	Real const& operator()(int i, int j) const { return coords_[i][j]; }

	XForm transpose() const;
	Real determinant() const;
	XForm inverse() const;

	template<int Dim2>
	XForm<Real, Dim2> cut() const;
private:
	Real coords_[Dim][Dim];
};

template<class Real, int Dim>
Point3D<Real> operator*(XForm<Real, Dim> const& f, Point3D<Real> const& p);

template<class Real, int Dim>
XForm<Real, Dim> operator*(XForm<Real, Dim> const& f1, XForm<Real, Dim> const& f2);

struct Triangle {
	double p[3][3];
};

struct CoredPointIndex {
	int index;
	char inCore;
};

struct TriangleIndex {
	std::array<int, 3> idx;
};

struct CoredVertexIndex {
	int idx;
	bool inCore;
};

class BufferedReadWriteFile {
public:
	BufferedReadWriteFile();
	~BufferedReadWriteFile();
	void reset();

	template<class T>
	bool write(T const& v) { return write(&v, sizeof(T)); }
	template<class T>
	bool read(T& v) { return read(&v, sizeof(T)); }

	template<class T>
	bool write(std::vector<T> const& v);
	template<class T>
	bool read(std::vector<T>& v);
private:
	bool write(void const* data, size_t size);
	bool read(void* data, size_t size);
private:
	FILE* fp_;
	char* buffer_;
	size_t buffer_index_;
	size_t buffer_size_;
};

template<class Vertex>
class CoredFileMeshData {
public:
	CoredFileMeshData(): out_of_core_points_file_(new BufferedReadWriteFile()),
		polygons_file_(new BufferedReadWriteFile()),
		out_of_core_points_count_(0), polygon_count_(0) { }
	~CoredFileMeshData() { delete out_of_core_points_file_; delete polygons_file_; }

	void resetIterator();

	void addInCorePoint(Vertex const& p) { return in_core_points_.push_back(p); }
	Vertex const& inCorePoints(int idx) { return in_core_points_[idx]; }
	int inCorePointCount() { return in_core_points_.size(); }

	int addOutOfCorePoint(Vertex const& p);
	bool nextOutOfCorePoint(Vertex& p) { return out_of_core_points_file_->read(p); }
	int outOfCorePointCount() { return out_of_core_points_count_; }

	int addPolygon(std::vector<CoredVertexIndex> const& vertices);
	bool nextPolygon(std::vector<CoredVertexIndex>& vs) { return polygons_file_->read(vs); }
	int polygonCount() { return polygon_count_; }
private:
	std::vector<Vertex> in_core_points_;
	// This fields are here to ensure this structure has the exact same memory layout as
	// before refactoring. Apparently it's used somewhere.
	// TODO: Find and eliminate.
	char pointFileName[1024];
	char polygonFileName[1024];
	BufferedReadWriteFile* out_of_core_points_file_;
	BufferedReadWriteFile* polygons_file_;
	int out_of_core_points_count_;
	int polygon_count_;
};

#include "Geometry.inl"
