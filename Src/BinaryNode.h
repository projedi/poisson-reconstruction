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

template<class Real>
class BinaryNode {
public:
	static int CenterCount(int depth) { return  1 << depth; }
	static int CumulativeCenterCount(int maxDepth) { return (1 << (maxDepth + 1)) - 1; }
	static int CenterIndex(int depth, int offset) { return (1 << depth) + offset - 1; }

	static int CornerCount(int depth) { return (1 << depth) + 1; }
	static int CumulativeCornerCount(int maxDepth)
		{ return (1 << (maxDepth + 1)) + maxDepth; }
	static int CornerIndex(int depth, int offset)
		{ return (1 << depth) + offset + depth; }

	static int CornerIndex(int maxDepth, int depth, int offset, int forwardCorner)
		{ return (offset + forwardCorner) << (maxDepth - depth); }
	static Real CornerIndexPosition(int index, int maxDepth)
		{ return Real(index) / (1 << maxDepth); }
	static Real Width(int depth) { return Real(1.0 / (1 << depth)); }

	static std::pair<Real, Real> CenterAndWidth(std::pair<int, int> const& dao);

	static std::pair<Real, Real> CenterAndWidth(int idx)
		{ return CenterAndWidth(DepthAndOffset(idx)); }

	static std::pair<int, int> DepthAndOffset(int idx);
};

template<class Real>
std::pair<Real, Real> BinaryNode<Real>::CenterAndWidth(
		std::pair<int, int> const& dao) {
	Real width = 1.0 / (1 << dao.first);
	return std::make_pair((0.5 + dao.second) * width, width);
}

template<class Real>
std::pair<int, int> BinaryNode<Real>::DepthAndOffset(int idx) {
	int i = idx + 1;
	std::pair<int, int> dao;
	// MSVC_2010_FIX - no idea what's going on.
	dao.first = 0;
	while(i) {
		i >>= 1;
		++dao.first;
	}
	// MSVC_2010_FIX - no idea what's going on.
	--dao.first;
	dao.second = (idx + 1) - (1 << dao.first);
	return dao;
}
