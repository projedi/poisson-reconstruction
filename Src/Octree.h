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

#include "Allocator.h"
#include "BinaryNode.h"
#include "MarchingCubes.h"

namespace octree_internals {

template<int N, class T>
class Neighbors {
public:
	T* neighbors[N][N][N];
	Neighbors(): neighbors{} { }
	void clear();
};

template<class OctNode, class Neighbors3, class Neighbors5>
class NeighborKey3 {
public:
	NeighborKey3(): neighbors_(nullptr), _depth(-1) { }
	NeighborKey3(NeighborKey3 const& nKey3);
	~NeighborKey3() { if(neighbors_) delete[] neighbors_; }

	Neighbors3* neighbors() const { return neighbors_; }

	void set(int depth);

	Neighbors3& setNeighbors(OctNode* node, bool flags[3][3][3] = (bool[3][3][3]){
		{ { true, true, true }, { true, true, true }, { true, true, true } },
		{ { true, true, true }, { true, true, true }, { true, true, true } },
		{ { true, true, true }, { true, true, true }, { true, true, true } }});

	Neighbors3& getNeighbors3(OctNode* node, int minDepth = 0);
	Neighbors5 getNeighbors5(OctNode* node);
private:
	Neighbors3& collectNeighbors(OctNode*, int minDepth, bool flags[3][3][3], bool doReset,
			std::function<void(OctNode*)> const& emptyChildrenCallback);
private:
	Neighbors3* neighbors_;
	int _depth;
};

}

#define DIMENSION 3

template<class NodeData, class Real>
class OctNode {
public:
	typedef octree_internals::Neighbors<3, OctNode> Neighbors3;
	typedef octree_internals::Neighbors<3, OctNode const> ConstNeighbors3;
	typedef octree_internals::Neighbors<5, OctNode> Neighbors5;
	typedef octree_internals::Neighbors<5, OctNode const> ConstNeighbors5;

	typedef octree_internals::NeighborKey3<OctNode, Neighbors3, Neighbors5> NeighborKey3;
	typedef octree_internals::NeighborKey3<OctNode const, ConstNeighbors3, ConstNeighbors5>
		ConstNeighborKey3;
public:
	static int const DepthShift = 5;
	static int const OffsetShift = (sizeof(long long) * 8 - DepthShift) / 3;
	static int const OffsetShift1 = DepthShift;
	static int const OffsetShift2 = OffsetShift1 + OffsetShift;
	static int const OffsetShift3 = OffsetShift2 + OffsetShift;
	static int const DepthMask = (1 << DepthShift) - 1;
	static int const OffsetMask = (1 << OffsetShift) - 1;;
	NodeData nodeData;
public:
	static void SetAllocator(int blockSize);

	template<class NodeAdjacencyFunction>
	static void ProcessFixedDepthNodeAdjacentNodes(int maxDepth,
			OctNode* node1, int width1, OctNode* node2, int width2,
			int depth, NodeAdjacencyFunction* F, int processCurrent = 1);

	static int CornerIndex(Point3D<Real> const& center, Point3D<Real> const& p);
public:
	OctNode(): parent_(nullptr), children_(nullptr), _depthAndOffset(0) { }
	~OctNode() { if(!UseAlloc && children_) delete[] children_; }

	OctNode* parent() const { return parent_; }
	OctNode*& children() const { return const_cast<OctNode*&>(children_); }

	bool initChildren();

	void depthAndOffset(int& depth, int offset[3]) const; 
	int depth() const { return _depthAndOffset & DepthMask; }
	void centerAndWidth(Point3D<Real>& center, Real& width) const;

	size_t leaves() const;
	size_t nodes() const;
	int maxDepth() const;

	OctNode* nextLeaf(OctNode* currentLeaf = nullptr);

	OctNode* nextNode(OctNode* currentNode = nullptr);

	OctNode* nextBranch(OctNode* current) const;

	void setFullDepth(int maxDepth);

	template<class NodeAdjacencyFunction>
	void processNodeFaces(OctNode const* node, NodeAdjacencyFunction& F, int fIndex,
			bool processCurrent = true) const;
private:
	template<class NodeAdjacencyFunction>
	static void ProcessFixedDepthNodeAdjacentNodes(int dx, int dy, int dz,
			OctNode* node1, int radius1, OctNode* node2, int radius2,
			int width2, int depth, NodeAdjacencyFunction* F, bool processCurrent = true);

	template<class NodeAdjacencyFunction>
	static void __ProcessFixedDepthNodeAdjacentNodes(int dx, int dy, int dz,
			OctNode* node1, int radius1, OctNode* node2, int radius2,
			int cWidth2, int depth, NodeAdjacencyFunction* F);

	// This is made private because the division by two has been pulled out.
	static bool Overlap(int c1, int c2, int c3, int dWidth);
	static int ChildOverlap(int dx, int dy, int dz, int d, int cRadius2);

	static unsigned long long Index(int depth, int const offset[3]);
private:
	template<class NodeAdjacencyFunction>
	void processNodeFaces(OctNode const* node, NodeAdjacencyFunction& F,
			std::tuple<int, int, int, int> cIndex) const;

	int width(int maxDepth) const { return 1 << (maxDepth - depth()); }
	void centerIndex(int maxDepth, int index[DIMENSION]) const;
private:
	static bool UseAlloc;
	static Allocator<OctNode> allocator;

	OctNode* parent_;
	OctNode* children_;
	unsigned long long _depthAndOffset;
};

#include "Octree.inl"
