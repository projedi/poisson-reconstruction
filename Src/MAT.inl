/*
Copyright (c) 2007, Michael Kazhdan
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

template <class Real>
void MinimalAreaTriangulation<Real>::GetTriangulation(
		std::vector<Point3D<Real>> const& vertices, std::vector<TriangleIndex>& triangles) {
	if(vertices.size() == 3) {
		triangles.clear();
		triangles.push_back({ { { 0, 1, 2 } } });
		return;
	}
	if(vertices.size() == 4) {
		TriangleIndex tIndex[2][2] {
			{ { { { 0, 1, 2 } } }, { { { 2, 3, 0 } } } },
			{ { { { 0, 1, 3 } } }, { { { 3, 1, 2 } } } }
		};
		Real area[2]{};

		for(int i = 0; i != 2; ++i)
			for(int j = 0; j != 2; ++j)
				area[i] += (Real)Length(CrossProduct(
					vertices[tIndex[i][j].idx[1]] - vertices[tIndex[i][j].idx[0]],
					vertices[tIndex[i][j].idx[2]] - vertices[tIndex[i][j].idx[0]]));
		triangles.clear();
		int i = area[0] > area[1] ? 1 : 0;
		triangles.assign(tIndex[i], tIndex[i] + 2);
		return;
	}
	data_.clear();
	data_.resize(vertices.size() * vertices.size(), { -1, -1 });
	GetArea(0, 1, vertices);
	triangles.clear();
	GetTriangulation(0, 1, vertices, triangles);
}

template <class Real>
Real MinimalAreaTriangulation<Real>::GetArea(std::vector<Point3D<Real>> const& vertices) {
	data_.clear();
	data_.resize(vertices.size() * vertices.size(), { -1, -1 });
	return GetArea(0, 1, vertices);
}

template<class Real>
void MinimalAreaTriangulation<Real>::GetTriangulation(int i, int j,
		std::vector<Point3D<Real>> const& vertices, std::vector<TriangleIndex>& triangles) {
	if((i < j && i + vertices.size() <= (size_t)j + 1) || i == j || i == j + 1) return;

	int m = data_[i * vertices.size() + j].midPoint;
	if(m >= 0) {
		triangles.push_back({ { { i, j, m } } });
		GetTriangulation(i, m, vertices, triangles);
		GetTriangulation(m, j, vertices, triangles);
	}
}

template<class Real>
Real MinimalAreaTriangulation<Real>::GetArea(size_t i, size_t j,
		std::vector<Point3D<Real>> const& vertices) {
	size_t eCount = vertices.size();
	size_t idx = i * eCount + j;
	size_t ii = i < j ? i + eCount : i;
	if(j + 1 >= ii) {
		data_[idx].bestTriangulation = 0;
		return 0;
	}
	if(data_[idx].midPoint != -1) return data_[idx].bestTriangulation;

	Real a = FLT_MAX;
	int mid = -1;
	for(size_t r = j + 1; r < ii; ++r) {
		size_t rr = r % eCount;
		Real temp = Real(Length(CrossProduct(vertices[i] - vertices[rr],
			vertices[j] - vertices[rr])));
		size_t idx1 = i * eCount + rr;
		size_t idx2 = rr * eCount + j;
		if(data_[idx1].bestTriangulation >= 0) {
			temp += data_[idx1].bestTriangulation;
			if(temp > a) continue;
			temp += data_[idx2].bestTriangulation > 0 ?
				data_[idx2].bestTriangulation : GetArea(rr, j, vertices);
		} else {
			temp += data_[idx2].bestTriangulation >= 0 ?
				data_[idx2].bestTriangulation : GetArea(rr, j, vertices);
			if(temp > a) continue;
			temp += GetArea(i, rr, vertices);
		}

		if(temp < a) {
			a = temp;
			mid = rr;
		}
	}
	data_[idx] = { a, mid };

	return a;
}
