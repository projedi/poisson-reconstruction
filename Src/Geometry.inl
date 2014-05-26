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

/////////////
// Point3D //
/////////////

template<class Real>
Point3D<Real>& operator+=(Point3D<Real>& p1, Point3D<Real> const& p2) {
	for(int i = 0; i != 3; ++i) p1[i] += p2[i];
	return p1;
}

template<class Real>
Point3D<Real>& operator*=(Point3D<Real>& p, Real r) {
	for(int i = 0; i != 3; ++i) p[i] *= r;
	return p;
}

template<class Real>
Point3D<Real>& operator/=(Point3D<Real>& p, Real r) {
	for(int i = 0; i != 3; ++i) p[i] /= r;
	return p;
}

template <class Real>
Point3D<Real> CrossProduct(Point3D<Real> const& p1, Point3D<Real> const& p2) {
	return Point3D<Real>(p1[1] * p2[2] - p1[2] * p2[1], -p1[0] * p2[2] + p1[2] * p2[0],
		p1[0] * p2[1] - p1[1] * p2[0]);
}

///////////
// XForm //
///////////

template<class Real, int Dim>
XForm<Real, Dim> XForm<Real, Dim>::Identity() {
	XForm xForm;
	for(int i = 0; i != Dim; ++i) xForm(i, i) = 1;
	return xForm;
}

template<class Real, int Dim>
XForm<Real, Dim> XForm<Real, Dim>::transpose() const {
	XForm xForm;
	for(int i = 0; i != Dim; ++i)
		for(int j = 0; j != Dim; ++j)
			xForm(i, j) = coords_[j][i];
	return xForm;
}

template<class Real>
Real subDeterminant(XForm<Real, 3> const& f, int i, int j) {
	int i1 = (i + 1) % 3;
	int i2 = (i + 2) % 3;
	int j1 = (j + 1) % 3;
	int j2 = (j + 2) % 3;
	Real m = (i + j) % 2 ? -1 : 1;
	return m * (f(i1, j1) * f(i2, j2) - f(i1, j2) * f(i2, j1));
}

template<class Real>
Real subDeterminant(XForm<Real, 4> const& f, int i, int j) {
	XForm<Real, 3> xForm;
	int ii[] { (i + 1) % 4, (i + 2) % 4, (i + 3) % 4 };
	int jj[] { (j + 1) % 4, (j + 2) % 4, (j + 3) % 4 };
	for(int _i = 0; _i != 3; ++_i)
		for(int _j = 0; _j != 3; ++_j)
			xForm(_i, _j) = f(ii[_i], jj[_j]);
	return xForm.determinant();
}

template<class Real, int Dim>
Real XForm<Real, Dim>::determinant() const {
	Real r = 0;
	for(int i = 0; i != Dim; ++i)
		r += i % 2 ? -coords_[i][0] * subDeterminant(*this, i, 0) :
			coords_[i][0] * subDeterminant(*this, i, 0);
	return r;
}

template<class Real, int Dim>
XForm<Real, Dim> XForm<Real, Dim>::inverse() const {
	XForm xForm;
	Real d = determinant();
	for(int i = 0; i != Dim; ++i)
		for(int j = 0; j != Dim; ++j)
			xForm.coords_[j][i] = (i + j) % 2 ? -subDeterminant(*this, i, j) / d :
				subDeterminant(*this, i, j) / d;
	return xForm;
}

template<class Real, int Dim>
Point3D<Real> operator*(XForm<Real, Dim> const& f, Point3D<Real> const& p) {
	Point3D<Real> q;
	for(int i = 0; i != 3; ++i) {
		int j;
		for(j = 0; j != 3; ++j)
			q[i] += f(j, i) * p[j];
		for(; j < Dim; ++j)
			q[i] += f(j, i);
	}
	return q;
}

template<class Real, int Dim>
XForm<Real, Dim> operator*(XForm<Real, Dim> const& f1, XForm<Real, Dim> const& f2) {
	XForm<Real, Dim> n;
	for(int i = 0; i != Dim; ++i)
		for(int j = 0; j != Dim; ++j)
			for(int k = 0; k != Dim; ++k)
				n(i, j) += f2(i, k) * f1(k, j);
	return n;
}

///////////////////////////
// BufferedReadWriteFile //
///////////////////////////

template<class T>
bool BufferedReadWriteFile::write(std::vector<T> const& vs) {
	if(!write((int)vs.size())) return false;
	return write(vs.data(), sizeof(T) * vs.size());
}

template<class T>
bool BufferedReadWriteFile::read(std::vector<T>& vs) {
	int pSize;
	if(!read(pSize)) return false;
	vs.resize(pSize);
	return read(vs.data(), sizeof(T) * vs.size());
}

///////////////////////
// CoredFileMeshData //
///////////////////////

template<class Vertex>
void CoredFileMeshData<Vertex>::resetIterator() {
	out_of_core_points_file_->reset();
	polygons_file_->reset();
}

template<class Vertex>
int CoredFileMeshData<Vertex>::addOutOfCorePoint(Vertex const& p) {
	int sz;
#pragma omp critical (CoredFileMeshData_addOutOfCorePoint)
	{
		out_of_core_points_file_->write(p);
		sz = out_of_core_points_count_++;
	}
	return sz;
}

template<class Vertex>
int CoredFileMeshData<Vertex>::addPolygon(std::vector<CoredVertexIndex> const& vertices) {
	int sz;
#pragma omp critical (CoredFileMeshData_addPolygon)
	{
		polygons_file_->write(vertices);
		sz = polygon_count_++;
	}
	return sz;
}
