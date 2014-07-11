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

#include <stdlib.h>
#include <math.h>
#include <algorithm>

/////////////
// OctNode //
/////////////

template<class NodeData, class Real>
bool OctNode<NodeData, Real>::UseAlloc = false;

template<class NodeData, class Real>
Allocator<OctNode<NodeData, Real> > OctNode<NodeData, Real>::allocator;

template<class NodeData, class Real>
void OctNode<NodeData, Real>::SetAllocator(int blockSize) {
	if(blockSize > 0) {
		UseAlloc = true;
		allocator.set(blockSize);
	} else UseAlloc = false;
}

template<class NodeData, class Real>
void OctNode<NodeData, Real>::setFullDepth(int maxDepth) {
	if(!maxDepth) return;
	if(!children_) initChildren();
	for(int i = 0; i != 8; ++i) children_[i].setFullDepth(maxDepth - 1);
}

template<class NodeData, class Real>
bool OctNode<NodeData, Real>::initChildren() {
	if(UseAlloc) children_ = allocator.newElements(Cube::CORNERS);
	else {
		if(children_) delete[] children_;
		children_ = new OctNode[Cube::CORNERS];
	}
	if(!children_) {
		std::cerr << "[ERROR] Failed to initialize children_ in OctNode::initChildren" <<
			std::endl;
		std::exit(1);
	}
	int d;
	int off[3];
	depthAndOffset(d, off);
	for(int i = 0; i != 2; ++i)
		for(int j = 0; j != 2; ++j)
			for(int k = 0; k != 2; ++k) {
				int idx = Cube::CornerIndex(i, j, k);
				children_[idx].parent_ = this;
				children_[idx].children_ = nullptr;
				int off2[3] = { (off[0] << 1) + i, (off[1] << 1) + j, (off[2] << 1) + k };
				children_[idx]._depthAndOffset = Index(d + 1, off2);
			}
	return true;
}

template<class NodeData, class Real>
unsigned long long OctNode<NodeData, Real>::Index(int depth, int const offset[3]) {
	unsigned long long idx = 0;
	idx |= ((unsigned long long)depth) & DepthMask;
	idx |= (((unsigned long long)offset[0]) & OffsetMask) << OffsetShift1;
	idx |= (((unsigned long long)offset[1]) & OffsetMask) << OffsetShift2;
	idx |= (((unsigned long long)offset[2]) & OffsetMask) << OffsetShift3;
	return idx;
}

template<class NodeData, class Real>
void OctNode<NodeData, Real>::depthAndOffset(int& depth, int offset[3]) const {
	depth = _depthAndOffset & DepthMask;
	offset[0] = (_depthAndOffset >> OffsetShift1) & OffsetMask;
	offset[1] = (_depthAndOffset >> OffsetShift2) & OffsetMask;
	offset[2] = (_depthAndOffset >> OffsetShift3) & OffsetMask;
}

template<class NodeData, class Real>
void OctNode<NodeData, Real>::centerAndWidth(Point3D<Real>& center, Real& width) const {
	int depth;
	int offset[3];
	depthAndOffset(depth, offset);
	width = 1.0 / (1 << depth);
	for(int dim = 0; dim != 3; ++dim)
		center.coords[dim] = (0.5 + offset[dim]) * width;
}

template<class NodeData, class Real>
int OctNode<NodeData, Real>::maxDepth() const {
	if(!children_) return 0;
	int c = 0;
	for(unsigned i = 0; i != Cube::CORNERS; ++i) {
		int d = children_[i].maxDepth();
		if(d > c) c = d;
	}
	return c + 1;
}

template<class NodeData, class Real>
size_t OctNode<NodeData, Real>::nodes() const {
	if(!children_) return 1;
	size_t c = 0;
	for(unsigned i = 0; i != Cube::CORNERS; ++i) c += children_[i].nodes();
	return c + 1;
}

template<class NodeData, class Real>
size_t OctNode<NodeData, Real>::leaves() const {
	if(!children_) return 1;
	size_t c = 0;
	for(unsigned i = 0; i != Cube::CORNERS; ++i) c += children_[i].leaves();
	return c;
}

template<class NodeData, class Real>
OctNode<NodeData, Real>* OctNode<NodeData, Real>::nextBranch(OctNode* current) const {
	if(!current->parent_ || current == this) return nullptr;
	if(current - current->parent_->children_ == Cube::CORNERS - 1)
		return nextBranch(current->parent_);
	return current + 1;
}

template<class NodeData, class Real>
OctNode<NodeData, Real>* OctNode<NodeData, Real>::nextLeaf(OctNode* current) {
	if(!current) {
		OctNode* res = this;
		while(res->children_) res = res->children_;
		return res;
	}
	if(current->children_) return current->nextLeaf();
	OctNode* res = nextBranch(current);
	return res ? res->nextLeaf() : nullptr;
}

template<class NodeData, class Real>
OctNode<NodeData, Real>* OctNode<NodeData, Real>::nextNode(OctNode* current) {
	if(!current) return this;
	if(current->children_) return current->children_;
	return nextBranch(current);
}

template<class NodeData, class Real>
template<class NodeAdjacencyFunction>
void OctNode<NodeData, Real>::processNodeFaces(OctNode const* node,
		NodeAdjacencyFunction const& F, int fIndex, bool processCurrent) const {
	if(processCurrent) F(this, node);
	int cs[4];
	Cube::FaceCorners(fIndex, cs[0], cs[1], cs[2], cs[3]);
	if(children_) processNodeFaces(node, F, cs);
}

template<class NodeData, class Real>
template<class NodeAdjacencyFunction>
void OctNode<NodeData, Real>::processNodeFaces(OctNode const* node,
		NodeAdjacencyFunction const& F, int cIndex[4]) const {
	F(&children_[cIndex[0]], node);
	F(&children_[cIndex[1]], node);
	F(&children_[cIndex[2]], node);
	F(&children_[cIndex[3]], node);

	if(children_[cIndex[0]].children_) children_[cIndex[0]].processNodeFaces(node, F, cIndex);
	if(children_[cIndex[1]].children_) children_[cIndex[1]].processNodeFaces(node, F, cIndex);
	if(children_[cIndex[2]].children_) children_[cIndex[2]].processNodeFaces(node, F, cIndex);
	if(children_[cIndex[3]].children_) children_[cIndex[3]].processNodeFaces(node, F, cIndex);
}

template<class NodeData, class Real>
template<class NodeAdjacencyFunction>
void OctNode<NodeData, Real>::ProcessFixedDepthNodeAdjacentNodes(int maxDepth,
		OctNode* node1, int width1, OctNode* node2, int width2,
		int depth, NodeAdjacencyFunction const& F, int processCurrent) {
	int c1[3];
	int c2[3];
	node1->centerIndex(maxDepth + 1, c1);
	node2->centerIndex(maxDepth + 1, c2);
	int w1 = node1->width(maxDepth + 1);
	int w2 = node2->width(maxDepth + 1);

	ProcessFixedDepthNodeAdjacentNodes(c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2],
			node1, (width1 * w1) >> 1, node2, (width2 * w2) >> 1, w2, depth, F,
			processCurrent);
}

template<class NodeData, class Real>
template<class NodeAdjacencyFunction>
void OctNode<NodeData, Real>::ProcessFixedDepthNodeAdjacentNodes(
		int dx, int dy, int dz, OctNode* node1, int radius1, OctNode* node2, int radius2,
		int width2, int depth, NodeAdjacencyFunction const& F, bool processCurrent) {
	int d = node2->depth();
	if(d > depth) return;
	if(!Overlap(dx, dy, dz, radius1 + radius2)) return;
	if(d == depth) {
		if(processCurrent) F(node2, node1);
	} else {
		if(!node2->children_) return;
		__ProcessFixedDepthNodeAdjacentNodes(-dx, -dy, -dz,
				node1, radius1, node2, radius2, width2 / 2, depth - 1, F);
	}
}

template<class NodeData, class Real>
template<class NodeAdjacencyFunction>
void OctNode<NodeData, Real>::__ProcessFixedDepthNodeAdjacentNodes(
		int dx, int dy, int dz, OctNode* node1, int radius1, OctNode* node2, int radius2,
		int cWidth2, int depth, NodeAdjacencyFunction const& F) {
	int cWidth = cWidth2 >> 1;
	int radius = radius2 >> 1;
	int o = ChildOverlap(dx, dy, dz, radius1 + radius, cWidth);
	if(!o) return;
	int dx1 = dx - cWidth;
	int dx2 = dx + cWidth;
	int dy1 = dy - cWidth;
	int dy2 = dy + cWidth;
	int dz1 = dz - cWidth;
	int dz2 = dz + cWidth;
	if(node2->depth() == depth) {
		if(o & 1) F(&node2->children_[0], node1);
		if(o & 2) F(&node2->children_[1], node1);
		if(o & 4) F(&node2->children_[2], node1);
		if(o & 8) F(&node2->children_[3], node1);
		if(o & 16) F(&node2->children_[4], node1);
		if(o & 32) F(&node2->children_[5], node1);
		if(o & 64) F(&node2->children_[6], node1);
		if(o & 128) F(&node2->children_[7], node1);
	} else {
		if(o & 1)
			if(node2->children_[0].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx1, dy1, dz1, node1, radius1,
						&node2->children_[0], radius, cWidth, depth, F);
		if(o & 2)
			if(node2->children_[1].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx2, dy1, dz1, node1, radius1,
						&node2->children_[1], radius, cWidth, depth, F);
		if(o & 4)
			if(node2->children_[2].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx1, dy2, dz1, node1, radius1,
						&node2->children_[2], radius, cWidth, depth, F);
		if(o & 8)
			if(node2->children_[3].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx2, dy2, dz1, node1, radius1,
						&node2->children_[3], radius, cWidth, depth, F);
		if(o & 16)
			if(node2->children_[4].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx1, dy1, dz2, node1, radius1,
						&node2->children_[4], radius, cWidth, depth, F);
		if(o & 32)
			if(node2->children_[5].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx2, dy1, dz2, node1, radius1,
						&node2->children_[5], radius, cWidth, depth, F);
		if(o & 64)
			if(node2->children_[6].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx1, dy2, dz2, node1, radius1,
						&node2->children_[6], radius, cWidth, depth, F);
		if(o & 128)
			if(node2->children_[7].children_)
				__ProcessFixedDepthNodeAdjacentNodes(dx2, dy2, dz2, node1, radius1,
						&node2->children_[7], radius, cWidth, depth, F);
	}
}

template<class NodeData, class Real>
int OctNode<NodeData, Real>::ChildOverlap(int dx, int dy, int dz, int d, int cRadius2) {
	int w1 = d - cRadius2;
	int w2 = d + cRadius2;
	int overlap = 0;

	int test = 0;
	int test1 = 0;
	if(dx < w2 && dx > -w1) test = 1;
	if(dx < w1 && dx > -w2) test |= 2;

	if(!test) return 0;
	if(dz < w2 && dz > -w1) test1 = test;
	if(dz < w1 && dz > -w2) test1 |= test << 4;

	if(!test1) return 0;
	if(dy < w2 && dy > -w1) overlap = test1;
	if(dy < w1 && dy > -w2) overlap |= test1 << 2;

	return overlap;
}

template<class NodeData, class Real>
int OctNode<NodeData, Real>::CornerIndex(Point3D<Real> const& center,
		Point3D<Real> const& p) {
	int cIndex = 0;
	if(p.coords[0] > center.coords[0]) cIndex |= 1;
	if(p.coords[1] > center.coords[1]) cIndex |= 2;
	if(p.coords[2] > center.coords[2]) cIndex |= 4;
	return cIndex;
}

template<class NodeData, class Real>
bool OctNode<NodeData, Real>::Overlap(int c1, int c2, int c3, int dWidth) {
	return !(c1 >= dWidth || c1 <= -dWidth || c2 >= dWidth || c2 <= -dWidth ||
			c3 >= dWidth || c3 <= -dWidth);
}

template<class NodeData, class Real>
void OctNode<NodeData, Real>::centerIndex(int maxDepth, int index[DIMENSION]) const {
	int d;
	int o[3];
	depthAndOffset(d, o);
	for(int i = 0; i != DIMENSION; ++i)
		index[i] = BinaryNode<Real>::CornerIndex(maxDepth, d + 1, o[i] << 1, 1);
}

namespace octree_internals {

template<int N, class T>
void Neighbors<N, T>::clear() {
	for(int i = 0; i != N; ++i)
		for(int j = 0; j != N; ++j)
			for(int k = 0; k != N; ++k)
				neighbors[i][j][k] = nullptr;
}

template<class OctNode>
template<class EmptyChildrenCallback>
typename NeighborKey3<OctNode>::Neighbors3& NeighborKey3<OctNode>::collectNeighbors(
		OctNode* node, int minDepth, bool flags[3][3][3], bool doReset,
		EmptyChildrenCallback const& emptyChildrenCallback) {
	int d = node->depth();
	if(d < minDepth) {
		std::cerr << "[ERROR] Node depth lower than min-depth: " << d <<
			" < " << minDepth << std::endl;
		std::exit(1);
	}
	if(doReset && node == neighbors_[d].neighbors[1][1][1]) {
		bool reset = false;
		for(int i = 0; i != 3; ++i)
			for(int j = 0; j != 3; ++j)
				for(int k = 0; k != 3; ++k)
					if(flags[i][j][k] && !neighbors_[d].neighbors[i][j][k]) reset = true;
		if(reset) neighbors_[d].neighbors[1][1][1] = nullptr;
	}
	if(node != neighbors_[d].neighbors[1][1][1]) {
		neighbors_[d].clear();

		if(d == minDepth)
			neighbors_[d].neighbors[1][1][1] = node;
		else {
			int x1;
			int y1;
			int z1;
			int idx = node->parent()->childIndex(node);
			Cube::FactorCornerIndex(idx, x1, y1, z1);
			int x2;
			int y2;
			int z2;
			Cube::FactorCornerIndex((~idx) & 7, x2, y2, z2);
			for(int i = 0; i != 2; ++i)
				for(int j = 0; j != 2; ++j)
					for(int k = 0; k != 2; ++k)
						neighbors_[d].neighbors[x2 + i][y2 + j][z2 + k] =
							node->parent()->child(Cube::CornerIndex(i, j, k));

			Neighbors3& temp = collectNeighbors(node->parent(), minDepth, flags,
					doReset, emptyChildrenCallback);

			// Set the neighbors from across the faces
			{
				int i = x1 << 1;
				if(temp.neighbors[i][1][1]) {
					if(flags[i][1][1] && !temp.neighbors[i][1][1]->hasChildren())
						emptyChildrenCallback(temp.neighbors[i][1][1]);
					if(temp.neighbors[i][1][1]->hasChildren())
						for(int j = 0; j != 2; ++j)
							for(int k = 0; k != 2; ++k)
								neighbors_[d].neighbors[i][y2 + j][z2 + k] =
									temp.neighbors[i][1][1]->child(Cube::CornerIndex(x2, j, k));
				}
			}
			{
				int j = y1 << 1;
				if(temp.neighbors[1][j][1]) {
					if(flags[1][j][1] && !temp.neighbors[1][j][1]->hasChildren())
						emptyChildrenCallback(temp.neighbors[1][j][1]);
					if(temp.neighbors[1][j][1]->hasChildren())
						for(int i = 0; i != 2; ++i)
							for(int k = 0; k != 2; ++k)
								neighbors_[d].neighbors[x2 + i][j][z2 + k] =
									temp.neighbors[1][j][1]->child(Cube::CornerIndex(i, y2, k));
				}
			}
			{
				int k = z1 << 1;
				if(temp.neighbors[1][1][k]) {
					if(flags[1][1][k] && !temp.neighbors[1][1][k]->hasChildren())
						emptyChildrenCallback(temp.neighbors[1][1][k]);
					if(temp.neighbors[1][1][k]->hasChildren())
						for(int i = 0; i != 2; ++i)
							for(int j = 0; j != 2; ++j)
								neighbors_[d].neighbors[x2 + i][y2 + j][k] =
									temp.neighbors[1][1][k]->child(Cube::CornerIndex(i, j, z2));
				}
			}

			// Set the neighbors from across the edges
			{
				int i = x1 << 1;
				int j = y1 << 1;
				if(temp.neighbors[i][j][1]) {
					if(flags[i][j][1] && !temp.neighbors[i][j][1]->hasChildren())
						emptyChildrenCallback(temp.neighbors[i][j][1]);
					if(temp.neighbors[i][j][1]->hasChildren())
						for(int k = 0; k != 2; ++k)
							neighbors_[d].neighbors[i][j][z2 + k] =
								temp.neighbors[i][j][1]->child(Cube::CornerIndex(x2, y2, k));
				}
			}
			{
				int i = x1 << 1;
				int k = z1 << 1;
				if(temp.neighbors[i][1][k]) {
					if(flags[i][1][k] && !temp.neighbors[i][1][k]->hasChildren())
						emptyChildrenCallback(temp.neighbors[i][1][k]);
					if(temp.neighbors[i][1][k]->hasChildren())
						for(int j = 0; j != 2; ++j)
							neighbors_[d].neighbors[i][y2 + j][k] =
								temp.neighbors[i][1][k]->child(Cube::CornerIndex(x2, j, z2));
				}
			}
			{
				int j = y1 << 1;
				int k = z1 << 1;
				if(temp.neighbors[1][j][k]) {
					if(flags[1][j][k] && !temp.neighbors[1][j][k]->hasChildren())
						emptyChildrenCallback(temp.neighbors[1][j][k]);
					if(temp.neighbors[1][j][k]->hasChildren())
						for(int i = 0; i != 2; ++i)
							neighbors_[d].neighbors[x2 + i][j][k] =
								temp.neighbors[1][j][k]->child(Cube::CornerIndex(i, y2, z2));
				}
			}

			// Set the neighbor from across the corner
			{
				int i = x1 << 1;
				int j = y1 << 1;
				int k = z1 << 1;
				if(temp.neighbors[i][j][k]) {
					if(flags[i][j][k] && !temp.neighbors[i][j][k]->hasChildren())
						emptyChildrenCallback(temp.neighbors[i][j][k]);
					if(temp.neighbors[i][j][k]->hasChildren())
						neighbors_[d].neighbors[i][j][k] =
							temp.neighbors[i][j][k]->child(Cube::CornerIndex(x2, y2, z2));
				}
			}
		}
	}
	return neighbors_[d];
}

template<class OctNode>
void setNeighborsFunction(OctNode* node) { node->initChildren(); }

// Note the assumption is that if you enable an edge, you also enable adjacent faces.
// And, if you enable a corner, you enable adjacent edges and faces.
template<class OctNode>
typename NeighborKey3<OctNode>::Neighbors3& NeighborKey3<OctNode>::setNeighbors(OctNode* node,
		bool flags[3][3][3]) {
	return collectNeighbors(node, 0, flags, true, setNeighborsFunction<OctNode>);
}

template<class OctNode>
void getNeighbors3Function(OctNode*) { return; }

template<class OctNode>
typename NeighborKey3<OctNode>::Neighbors3& NeighborKey3<OctNode>::getNeighbors3(OctNode* node,
		int minDepth) {
	bool flags[3][3][3] = { };
	return collectNeighbors(node, minDepth, flags, false, getNeighbors3Function<OctNode>);
}

struct RangeData {
	RangeData(int s, int e, int o): start(s), end(e), neighborOffset(o) { }
	int start;
	int end;
	int neighborOffset;
};

RangeData getRange(int c, int i) {
	switch(i) {
		case 0: return c % 2 ? RangeData(1, 2, -1) : RangeData(0, 2, 0);
		case 1: return c % 2 ? RangeData(0, 2, 1) : RangeData(0, 2, 2);
		case 2: return c % 2 ? RangeData(0, 2, 3) : RangeData(0, 1, 4);
		default:
			std::cerr << "[ERROR]: Invalid value of i in getRange lambda" << std::endl;
			exit(1);
	}
};

template<class OctNode>
typename NeighborKey3<OctNode>::Neighbors5 NeighborKey3<OctNode>::getNeighbors5(OctNode* node) {
	Neighbors5 neighbors;
	if(!node) return neighbors;
	if(!node->parent()) {
		neighbors.neighbors[2][2][2] = node;
		return neighbors;
	}
	OctNode const** nodeIter = reinterpret_cast<OctNode const**>(
			getNeighbors3(node->parent()).neighbors);
	int idx = node->parent()->childIndex(node);
	for(int i = 0; i != 3; ++i)
		for(int j = 0; j != 3; ++j)
			for(int k = 0; k != 3; ++k, ++nodeIter) {
				RangeData ri = getRange(idx, i);
				RangeData rj = getRange(idx / 2, j);
				RangeData rk = getRange(idx / 4, k);
				if(*nodeIter && (*nodeIter)->hasChildren())
					for(int ii = ri.start; ii != ri.end; ++ii)
						for(int jj = rj.start; jj != rj.end; ++jj)
							for(int kk = rk.start; kk != rk.end; ++kk)
								neighbors.neighbors[ii + ri.neighborOffset]
									[jj + rj.neighborOffset][kk + rk.neighborOffset] =
										(*nodeIter)->child(Cube::CornerIndex(ii, jj, kk));
			}
	return neighbors;
}

}
