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

#include "Geometry.h"

struct Square {
	static unsigned const CORNERS = 4;
	static int CornerIndex(int x, int y) { return (y << 1) | x; }
	static int AntipodalCornerIndex(int idx);
	static void FactorCornerIndex(int idx, int& x, int& y);
};

struct Cube {
	static unsigned const CORNERS = 8;
	static unsigned const EDGES = 12;
	static unsigned const NEIGHBORS = 6;

	static int CornerIndex(int x, int y, int z) { return (z << 2) | (y << 1) | x; }
	static int EdgeIndex(int orient, int i, int j) { return (orient << 2) | (j << 1) | i; }

	static void FactorCornerIndex(int idx, int& x, int& y, int& z);
	static void FactorEdgeIndex(int idx, int& orient, int& i, int& j);
	static void FactorFaceIndex(int idx, int& dir, int& offset);
	static void FactorFaceIndexXYZ(int idx, int& x, int& y, int& z);

	static int AntipodalCornerIndex(int idx);
	static int FaceReflectEdgeIndex(int idx, int faceIndex);
	static int FaceReflectFaceIndex(int idx, int faceIndex);
	static int EdgeReflectEdgeIndex(int edgeIndex);

	static int FaceAdjacentToEdges(int eIndex1, int eIndex2);
	static void FacesAdjacentToEdge(int eIndex, int& f1, int& f2);

	static void EdgeCorners(int idx, int& c1, int& c2);
	static void FaceCorners(int idx, int& c1, int& c2, int& c3, int& c4);
private:
	static int FaceIndex(int x, int y, int z);
};

class MarchingCubes {
public:
	static unsigned const MAX_TRIANGLES = 5;
	static int const edgeMask[1 << Cube::CORNERS];
	static int const triangles[1 << Cube::CORNERS][3 * MAX_TRIANGLES + 1];
	static int const cornerMap[Cube::CORNERS];

	static int AddTriangleIndices(int idx, int* isoIndices);

	template<class Real>
	static int GetIndex(Real const values[Cube::CORNERS], Real iso);

	static bool HasRoots(int mcIndex) { return mcIndex != 0 && mcIndex != 255; }
	static int HasEdgeRoots(int mcIndex, int edgeIndex);
};

template<class Real>
int MarchingCubes::GetIndex(Real const v[Cube::CORNERS], Real iso) {
	int idx = 0;
	if(v[Cube::CornerIndex(0, 0, 0)] < iso) idx |= 1;
	if(v[Cube::CornerIndex(1, 0, 0)] < iso) idx |= 2;
	if(v[Cube::CornerIndex(1, 1, 0)] < iso) idx |= 4;
	if(v[Cube::CornerIndex(0, 1, 0)] < iso) idx |= 8;
	if(v[Cube::CornerIndex(0, 0, 1)] < iso) idx |= 16;
	if(v[Cube::CornerIndex(1, 0, 1)] < iso) idx |= 32;
	if(v[Cube::CornerIndex(1, 1, 1)] < iso) idx |= 64;
	if(v[Cube::CornerIndex(0, 1, 1)] < iso) idx |= 128;
	return idx;
}
