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

#include <cassert>

#include "DumpOutput.h"
#include "Octree.h"
#include "time.h"
#include "MemoryUsage.h"
#include "PointStream.h"
#include "MAT.h"
#include "Util.h"

double const ITERATION_POWER = 1.0 / 3;
Real const MATRIX_ENTRY_EPSILON = 0;
Real const EPSILON = 1e-6;
Real const ROUND_EPS = 1e-5;

//////////////////
// TreeNodeData //
//////////////////

template<bool StoreDensity>
TreeNodeData<StoreDensity>::TreeNodeData():
	nodeIndex(-1),
	normalIndex(-1),
	constraint(0),
	solution(0),
	pointIndex(-1) {
	centerWeightContribution[0] = 0;
	if(StoreDensity)
		centerWeightContribution[1] = 0;
}

/////////////////////
// SortedTreeNodes //
/////////////////////

template<bool OutputDensity>
void SortedTreeNodes<OutputDensity>::set(TreeOctNode& root) {
	treeNodes.clear();
	maxDepth = root.maxDepth() + 1;
	nodeCount.resize(maxDepth + 1);
	treeNodes.resize(root.nodes());

	int startDepth = 0;
	nodeCount[0] = 0;
	nodeCount[1] = 1;
	treeNodes[0] = &root;
	for(TreeOctNode* node = root.nextNode(); node; node = root.nextNode(node))
		node->nodeData.nodeIndex = -1;
	for(int d = startDepth + 1; d != maxDepth; ++d) {
		nodeCount[d + 1] = nodeCount[d];
		for(int i = nodeCount[d - 1]; i != nodeCount[d]; ++i) {
			TreeOctNode* temp = treeNodes[i];
			if(temp->hasChildren())
				for(int c = 0; c != 8; ++c)
					treeNodes[nodeCount[d + 1]++] = temp->child(c);
		}
	}
	for(int i = 0; i != nodeCount[maxDepth]; ++i) treeNodes[i]->nodeData.nodeIndex = i;
}

template<bool OutputDensity>
void SortedTreeNodes<OutputDensity>::setCornerTable(CornerTableData& cData, TreeOctNode const* rootNode,
		int maxDepth, int threads) const {
	if(threads <= 0) threads = 1;
	cData.resizeOffsets(this->maxDepth, -1);
	// The vector of per-depth node spans
	std::vector<std::pair<int, int> > spans(this->maxDepth, std::pair<int, int>(-1, -1));
	int minDepth;
	int off[3];
	int start;
	int end = 0;
	if(rootNode) {
		rootNode->depthAndOffset(minDepth, off);
		start = end = rootNode->nodeData.nodeIndex;
	} else {
		start = 0;
		for(minDepth = 0; minDepth <= this->maxDepth; ++minDepth)
			if(nodeCount[minDepth + 1]){
				end = nodeCount[minDepth + 1] - 1;
				break;
			}
	}
	int nodeCount = 0;
	for(int d = minDepth; d <= maxDepth; ++d) {
		spans[d] = std::pair<int, int>(start, end + 1);
		cData.offsets(d) = nodeCount - spans[d].first;
		nodeCount += spans[d].second - spans[d].first;
		if(d < maxDepth) {
			while(start < end && !treeNodes[start]->hasChildren()) ++start;
			while(end > start && !treeNodes[end]->hasChildren()) --end;
			if(start == end && !treeNodes[start]->hasChildren()) break;
			start = treeNodes[start]->child(0)->nodeData.nodeIndex;
			end = treeNodes[end]->child(7)->nodeData.nodeIndex;
		}
	}

	cData.resizeTable(nodeCount);
	int count = 0;
	TreeConstNeighborKey3 neighborKey(maxDepth);
	std::vector<int> cIndices(nodeCount * Cube::CORNERS, 0);
	for(int d = minDepth; d <= maxDepth; ++d) {
#pragma omp parallel for num_threads(threads) firstprivate(neighborKey)
		for(int i = spans[d].first; i < spans[d].second; ++i) {
			TreeOctNode* node = treeNodes[i];
			if(d < maxDepth && node->hasChildren()) continue;
			typename TreeOctNode::ConstNeighbors3 const& neighbors =
				neighborKey.getNeighbors3(node, minDepth);
			for(unsigned c = 0; c != Cube::CORNERS; ++c)	{ // Iterate over the cell's corners
				bool cornerOwner = true;
				int x;
				int y;
				int z;
				Cube::FactorCornerIndex(c, x, y, z);
				unsigned ac = Cube::AntipodalCornerIndex(c); // The index of the node relative to the corner
				for(unsigned cc = 0; cc != Cube::CORNERS; ++cc) { // Iterate over the corner's cells
					int xx;
					int yy;
					int zz;
					Cube::FactorCornerIndex(cc, xx, yy, zz);
					xx += x;
					yy += y;
					zz += z;
					if(neighbors.at(xx, yy, zz) &&
							neighbors.at(xx, yy, zz)->nodeData.nodeIndex != -1)
						if(cc < ac || (d < maxDepth && neighbors.at(xx, yy, zz)->hasChildren())) {
							int _d;
							int _off[3];
							neighbors.at(xx, yy, zz)->depthAndOffset(_d, _off);
							_off[0] >>= d - minDepth;
							_off[1] >>= d - minDepth;
							_off[2] >>= d - minDepth;
							if(!rootNode || (_off[0] == off[0] && _off[1] == off[1] && _off[2] == off[2])) {
								cornerOwner = false;
								break;
							}
							else std::cerr << "[WARNING] How did we leave the subtree?" << std::endl;
						}
				}
				if(cornerOwner) {
					int myCount = (treeNodes[i]->nodeData.nodeIndex + cData.offsets(d)) * Cube::CORNERS + c;
					cIndices[myCount] = 1;
					TreeOctNode const* n = node;
					int d = n->depth();
					do {
						typename TreeOctNode::ConstNeighbors3 const& neighbors = neighborKey.neighbors(d);
						// Set all the corner indices at the current depth
						for(unsigned cc = 0; cc != Cube::CORNERS; ++cc) {
							int xx;
							int yy;
							int zz;
							Cube::FactorCornerIndex(cc, xx, yy, zz);
							xx += x;
							yy += y;
							zz += z;
							if(neighbors.at(xx, yy, zz) && neighbors.at(xx, yy, zz)->nodeData.nodeIndex != -1)
								cData[neighbors.at(xx, yy, zz)][Cube::AntipodalCornerIndex(cc)] = myCount;
						}
						// If we are not at the root and the parent() also has the corner
						if(d == minDepth || n != n->parent()->child(c)) break;
						n = n->parent();
						--d;
					} while(1);
				}
			}
		}
	}
	for(size_t i = 0; i != cIndices.size(); ++i)
		if(cIndices[i]) cIndices[i] = count++;
	for(int d = minDepth; d <= maxDepth; ++d)
#pragma omp parallel for num_threads(threads)
		for(int i = spans[d].first; i < spans[d].second; ++i)
			for(unsigned j = 0; j != Cube::CORNERS; ++j)
				cData[treeNodes[i]][j] = cIndices[cData[treeNodes[i]][j]];
	cData.setCount(count);
}

template<bool OutputDensity>
int SortedTreeNodes<OutputDensity>::getMaxCornerCount(int depth, int maxDepth, int threads) const {
	if(threads <= 0) threads = 1;
	int res = 1 << depth;

	std::vector<int> cornerCount(res * res * res, 0);
	TreeConstNeighborKey3 neighborKey(maxDepth);
#pragma omp parallel for num_threads(threads) firstprivate(neighborKey)
	for(int i = nodeCount[depth]; i < nodeCount[maxDepth + 1]; ++i) {
		TreeOctNode* node = treeNodes[i];
		int d;
		int off[3];
		node->depthAndOffset(d, off);
		if(d < maxDepth && node->hasChildren()) continue;

		typename TreeOctNode::ConstNeighbors3 const& neighbors = neighborKey.getNeighbors3(node, depth);
		for(unsigned c = 0; c != Cube::CORNERS; ++c) { // Iterate over the cell's corners
			bool cornerOwner = true;
			int x;
			int y;
			int z;
			unsigned ac = Cube::AntipodalCornerIndex(c); // The index of the node relative to the corner
			Cube::FactorCornerIndex(c, x, y, z);
			for(unsigned cc = 0; cc != Cube::CORNERS; ++cc) { // Iterate over the corner's cells
				int xx, yy, zz;
				Cube::FactorCornerIndex(cc, xx, yy, zz);
				xx += x;
				yy += y;
				zz += z;
				if(neighbors.at(xx, yy, zz) &&
						neighbors.at(xx, yy, zz)->nodeData.nodeIndex != -1)
					if(cc < ac || (d < maxDepth && neighbors.at(xx, yy, zz)->hasChildren())) {
						cornerOwner = false;
						break;
					}
			}
			if(cornerOwner)
#pragma omp atomic
				++cornerCount[((off[0] >> (d - depth)) * res * res) + ((off[1] >> (d - depth)) * res) +
					(off[2] >> (d - depth))];
		}
	}
	int maxCount = 0;
	for(int i = 0; i != res * res * res; ++i) maxCount = std::max(maxCount, cornerCount[i]);
	return maxCount;
}

template<bool OutputDensity>
void SortedTreeNodes<OutputDensity>::setEdgeTable(EdgeTableData& eData, TreeOctNode const* rootNode,
		int maxDepth, int threads) {
	if(threads <= 0) threads = 1;
	std::vector<std::pair<int, int> > spans(this->maxDepth, std::pair<int, int>(-1, -1));

	int minDepth;
	eData.resizeOffsets(this->maxDepth, -1);
	int start;
	int end = 0;
	if(rootNode) {
		minDepth = rootNode->depth();
		start = end = rootNode->nodeData.nodeIndex;
	} else {
		start = 0;
		for(minDepth = 0; minDepth <= this->maxDepth; ++minDepth)
			if(nodeCount[minDepth+1]) {
				end = nodeCount[minDepth + 1] - 1;
				break;
			}
	}

	int nodeCount = 0;
	for(int d = minDepth; d <= maxDepth; ++d) {
		spans[d] = std::pair<int, int>(start, end + 1);
		eData.offsets(d) = nodeCount - spans[d].first;
		nodeCount += spans[d].second - spans[d].first;
		if(d < maxDepth) {
			while(start < end && !treeNodes[start]->hasChildren()) ++start;
			while(end > start && !treeNodes[end]->hasChildren()) --end;
			if(start == end && !treeNodes[start]->hasChildren()) break;
			start = treeNodes[start]->child(0)->nodeData.nodeIndex;
			end = treeNodes[end]->child(7)->nodeData.nodeIndex;
		}
	}
	eData.resizeTable(nodeCount);
	std::vector<int> eIndices(nodeCount * Cube::EDGES, 0);
	int count = 0;
	TreeConstNeighborKey3 neighborKey(maxDepth);
	for(int d = minDepth; d <= maxDepth; ++d) {
#pragma omp parallel for num_threads(threads) firstprivate(neighborKey)
		for(int i = spans[d].first; i < spans[d].second; ++i) {
			TreeOctNode* node = treeNodes[i];
			typename TreeOctNode::ConstNeighbors3 const& neighbors = neighborKey.getNeighbors3(node, minDepth);

			for(unsigned e = 0; e != Cube::EDGES; ++e) {
				bool edgeOwner = true;
				int o;
				int _i;
				int _j;
				Cube::FactorEdgeIndex(e, o, _i, _j);
				unsigned ac = Square::AntipodalCornerIndex(Square::CornerIndex(_i, _j));
				for(unsigned cc = 0; cc != Square::CORNERS; ++cc) {
					int ii;
					int jj;
					int x = 0;
					int y = 0;
					int z = 0;
					Square::FactorCornerIndex(cc, ii, jj);
					ii += _i;
					jj += _j;
					switch(o) {
						case 0: x = 1; y = ii; z = jj; break;
						case 1: x = ii; y = 1; z = jj; break;
						case 2: x = ii; y = jj; z = 1; break;
					}
					if(neighbors.at(x, y, z) &&
							neighbors.at(x, y, z)->nodeData.nodeIndex != -1 && cc < ac) {
						edgeOwner = false;
						break;
					} 
				}
				if(edgeOwner) {
					int myCount = (treeNodes[i]->nodeData.nodeIndex + eData.offsets(d)) * Cube::EDGES + e;
					eIndices[myCount] = 1;
					// Set all edge indices
					for(unsigned cc = 0; cc != Square::CORNERS; ++cc) {
						int ii;
						int jj;
						int aii;
						int ajj;
						int x = 0;
						int y = 0;
						int z = 0;
						Square::FactorCornerIndex(cc, ii, jj);
						Square::FactorCornerIndex(Square::AntipodalCornerIndex(cc), aii, ajj);
						ii += _i;
						jj += _j;
						switch(o) {
							case 0: x = 1; y = ii; z = jj; break;
							case 1: x = ii; y = 1; z = jj; break;
							case 2: x = ii; y = jj; z = 1; break;
						}
						if(neighbors.at(x, y, z) && neighbors.at(x, y, z)->nodeData.nodeIndex != -1)
							eData[neighbors.at(x, y, z)][Cube::EdgeIndex(o, aii, ajj)] = myCount;
					}
				}
			}
		}
	}
	for(size_t i = 0; i != eIndices.size(); ++i) if(eIndices[i]) eIndices[i] = count++;
	for(int d = minDepth; d <= maxDepth; ++d)
#pragma omp parallel for num_threads(threads)
		for(int i = spans[d].first; i < spans[d].second; ++i)
			for(unsigned j = 0; j != Cube::EDGES; ++j)
				eData[treeNodes[i]][j] = eIndices[eData[treeNodes[i]][j]];
	eData.setCount(count);
}

template<bool OutputDensity>
int SortedTreeNodes<OutputDensity>::getMaxEdgeCount(TreeOctNode const*, int depth, int threads) const {
	if(threads <= 0) threads = 1;
	int res = 1 << depth;
	std::vector<int> edgeCount(res * res * res, 0);
	TreeConstNeighborKey3 neighborKey(maxDepth -1);
#pragma omp parallel for num_threads(threads) firstprivate(neighborKey)
	for(int i = nodeCount[depth]; i < nodeCount[maxDepth]; ++i) {
		TreeOctNode* node = treeNodes[i];
		typename TreeOctNode::ConstNeighbors3 const& neighbors = neighborKey.getNeighbors3(node, depth);
		int d;
		int off[3];
		node->depthAndOffset(d, off);

		for(unsigned e = 0; e != Cube::EDGES; ++e) {
			bool edgeOwner = true;
			int o;
			int i;
			int j;
			Cube::FactorEdgeIndex(e, o, i, j);
			unsigned ac = Square::AntipodalCornerIndex(Square::CornerIndex(i, j));
			for(unsigned cc = 0; cc != Square::CORNERS; ++cc) {
				int ii;
				int jj;
				int x = 0;
				int y = 0;
				int z = 0;
				Square::FactorCornerIndex(cc, ii, jj);
				ii += i;
				jj += j;
				switch(o) {
					case 0: x = 1; y = ii; z = jj; break;
					case 1: x = ii; y = 1; z = jj; break;
					case 2: x = ii; y = jj; z = 1; break;
				}
				if(neighbors.at(x, y, z) && neighbors.at(x, y, z)->nodeData.nodeIndex!=-1 && cc<ac) {
					edgeOwner = false;
					break;
				}
			}
			if(edgeOwner)
#pragma omp atomic
				++edgeCount[((off[0] >> (d - depth)) * res * res) + ((off[1] >> (d - depth)) * res) +
					(off[2] >> (d - depth))];
		}
	}
	int maxCount = 0;
	for(int i = 0; i != res * res * res; ++i) maxCount = std::max(maxCount, edgeCount[i]);
	return maxCount;
}

////////////
// Octree //
////////////
template<int Degree, bool OutputDensity>
size_t Octree<Degree, OutputDensity>::maxMemoryUsage_ = 0;

template<int Degree, bool OutputDensity>
double Octree<Degree, OutputDensity>::MemoryUsage() {
	size_t mem = MemoryInfo::Usage();
	if(mem > maxMemoryUsage_) maxMemoryUsage_ = mem;
	return mem;
}

template<int Degree, bool OutputDensity>
Octree<Degree, OutputDensity>::Octree(int threads, int maxDepth, BoundaryType boundaryType):
	threads_(threads),
	boundaryType_(boundaryType),
	radius_(0.5 + 0.5 * Degree),
	width_((int)((double)(radius_ + 0.5 - EPSILON) * 2)),
	constrainValues_(false) {
	if(boundaryType_ == BoundaryTypeNone) ++maxDepth;
	postDerivativeSmooth_ = (Real)1.0 / (1 << maxDepth);
	fData_.set(maxDepth, (BoundaryType)boundaryType);
}

template<int Degree, bool OutputDensity>
bool Octree<Degree, OutputDensity>::IsInset(TreeOctNode const* node) {
	int d;
	int off[3];
	node->depthAndOffset(d, off);
	int res = 1 << d;
	int o = 1 << (d - 2);
	return off[0] >= o && off[0] < res - o &&
		off[1] >= o && off[1] < res - o &&
		off[2] >= o && off[2] < res - o;
}

template<class TreeOctNode>
bool IsInsetSupported(TreeOctNode const* node) {
	int d;
	int off[3];
	node->depthAndOffset(d, off);
	int res = 1 << d;
	int o = 1 << (d - 2);
	return off[0] >= o && off[0] <= res - o &&
		off[1] >= o && off[1] <= res - o &&
		off[2] >= o && off[2] <= res - o;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::SplatOrientedPoint(TreeOctNode* node, Point3D<Real> const& position,
		Point3D<Real> const& normal, TreeNeighborKey3& neighborKey) {
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center, w);
	double width = w;
	double dx[DIMENSION][SPLAT_ORDER + 1];
	int off[3];
	for(int i = 0; i != 3; ++i) {
#if SPLAT_ORDER == 2
		off[i] = 0;
		double x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.5 * x + 0.5 * x * x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.75 - x * x;
		dx[i][2] = 1 - dx[i][1] - dx[i][0];
#elif SPLAT_ORDER == 1
		double x = (position[i] - center[i]) / width;
		off[i] = x < 0 ? 0 : 1;
		dx[i][0] = x < 0 ? -x : 1 - x;
		dx[i][1] = 1 - dx[i][0];
#else
		// There was also SPLAT_ORDER == 0 case but it fails on the original code and so is removed.
#error Splat order not supported
#endif
	}
	TreeNeighbors3& neighbors = neighborKey.setNeighbors(node);
	for(int i = off[0]; i <= off[0] + SPLAT_ORDER; ++i) {
		for(int j = off[1]; j <= off[1] + SPLAT_ORDER; ++j) {
			for(int k = off[2]; k <= off[2] + SPLAT_ORDER; ++k) {
				TreeOctNode* nnode = neighbors.at(i, j, k);
				if(nnode) {
					int idx = nnode->nodeData.normalIndex;
					if(idx < 0) {
						nnode->nodeData.nodeIndex = 0;
						idx = nnode->nodeData.normalIndex = normals_.size();
						normals_.push_back(Point3D<Real>());
					}
					normals_[idx] += normal * (Real)(dx[0][i] * dx[1][j] * dx[2][k]);
				}
			}
		}
	}
}

template<class T1, class T2>
constexpr T1 clamp(T1 const& x, T2 const& minVal, T2 const& maxVal) {
	return x < minVal ? minVal : x > maxVal ? maxVal : x;
}

template<int Degree, bool OutputDensity>
Real Octree<Degree, OutputDensity>::SplatOrientedPoint(Point3D<Real> const& position,
		Point3D<Real> const& normal, TreeNeighborKey3& neighborKey, int splatDepth, Real samplesPerNode,
		int minDepth, int maxDepth) {
	TreeOctNode* temp = &tree_;
	Point3D<Real> myCenter(0.5, 0.5, 0.5);
	Real myWidth = 1;

	while(temp->depth() < splatDepth) {
		if(!temp->hasChildren()) {
			std::cerr << "Octree<Degree>::SplatOrientedPoint error" << std::endl;
			return -1;
		}
		int cIndex = TreeOctNode::CornerIndex(myCenter, position);
		temp = temp->child(cIndex);
		myWidth /= 2;
		myCenter[0] += cIndex & 1 ? myWidth / 2 : -myWidth / 2;
		myCenter[1] += cIndex & 2 ? myWidth / 2 : -myWidth / 2;
		myCenter[2] += cIndex & 4 ? myWidth / 2 : -myWidth / 2;
	}
	Real weight;
	Real depth;
	GetSampleDepthAndWeight(temp, position, SplatOrientedPointGetNeighborsFunction(neighborKey),
			samplesPerNode, depth, weight);

	depth = clamp(depth, minDepth, maxDepth);
	int topDepth = clamp(lrint(std::ceil(depth)), minDepth, maxDepth);

	while(temp->depth() > topDepth) temp = temp->parent();

	while(temp->depth() < topDepth) {
		if(!temp->hasChildren()) temp->initChildren();
		int cIndex = TreeOctNode::CornerIndex(myCenter, position);
		temp = temp->child(cIndex);
		myWidth /= 2;
		myCenter[0] += cIndex & 1 ? myWidth / 2 : -myWidth / 2;
		myCenter[1] += cIndex & 2 ? myWidth / 2 : -myWidth / 2;
		myCenter[2] += cIndex & 4 ? myWidth / 2 : -myWidth / 2;
	}
	Real dx = 1 - (topDepth - depth);
	double width = 1.0 / (1 << temp->depth());
	Point3D<Real> n = normal * weight / (Real)std::pow(width, 3) * dx;
	SplatOrientedPoint(temp, position, n, neighborKey);
	if(std::abs(1 - dx) > EPSILON) {
		dx = 1 - dx;
		temp = temp->parent();
		width = 1.0 / (1 << temp->depth());
		n = normal * weight / (Real)std::pow(width, 3) * dx;
		SplatOrientedPoint(temp, position, n, neighborKey);
	}
	return weight;
}

template<int Degree, bool OutputDensity>
template<class OctNodeS, class GetNeighbors>
void Octree<Degree, OutputDensity>::GetSampleDepthAndWeight(OctNodeS* node, Point3D<Real> const& position,
		GetNeighbors const& getNeighbors, Real samplesPerNode, Real& depth, Real& weight) const {
	OctNodeS* temp = node;
	weight = (Real)1.0 / GetSampleWeight(temp, position, getNeighbors(temp));
	if(weight >= samplesPerNode)
		depth = temp->depth() +
				std::log((double)(weight / samplesPerNode)) / std::log((double)(1 << (DIMENSION - 1)));
	else {
		Real oldWeight;
		Real newWeight;
		oldWeight = newWeight = weight;
		while(newWeight < samplesPerNode && temp->parent()) {
			temp = temp->parent();
			oldWeight = newWeight;
			newWeight = (Real)1.0 / GetSampleWeight(temp, position, getNeighbors(temp));
		}
		depth = temp->depth() +
			std::log((double)(newWeight / samplesPerNode)) / std::log((double)(newWeight / oldWeight));
	}
	weight = std::pow((double)(1 << (DIMENSION - 1)), -depth);
}

template<int Degree, bool OutputDensity>
Real Octree<Degree, OutputDensity>::GetSampleWeight(TreeOctNode const* node, Point3D<Real> const& position,
		TreeConstNeighbors3& neighbors) const {
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center, w);
	double width = w;
	double dx[DIMENSION][3];
	for(int i = 0; i != DIMENSION; ++i) {
		double x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.5 * x + 0.5 * x * x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.75 - x * x;
		dx[i][2] = 1 - dx[i][1] - dx[i][0];
	}

	Real weight = 0;
	for(int i = 0; i != 3; ++i) {
		for(int j = 0; j != 3; ++j) {
			for(int k = 0; k != 3; ++k) {
				if(neighbors.at(i, j, k))
					weight += (Real)(dx[0][i] * dx[1][j] * dx[2][k] *
						neighbors.at(i, j, k)->nodeData.centerWeightContribution[0]);
			}
		}
	}
	return 1.0 / weight;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::UpdateWeightContribution(TreeOctNode* node,
		Point3D<Real> const& position, TreeNeighborKey3& neighborKey, Real weight) const {
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center, w);
	double width = w;
	double const SAMPLE_SCALE = 1 / (0.125 * 0.125 + 0.75 * 0.75 + 0.125 * 0.125);

	double dx[DIMENSION][3];
	for(int i = 0; i != DIMENSION; ++i) {
		double x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.5 * x + 0.5 * x * x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.75 - x * x;
		dx[i][2] = 1 - dx[i][1] - dx[i][0];
		// Note that we are splatting along a co-dimension one manifold, so uniform point samples
		// do not generate a unit sample weight.
		dx[i][0] *= SAMPLE_SCALE;
	}

	TreeNeighbors3& neighbors = neighborKey.setNeighbors(node);
	for(int i = 0; i != 3; ++i) {
		for(int j = 0; j != 3; ++j) {
			for(int k = 0; k != 3; ++k) {
				if(neighbors.at(i, j, k))
					neighbors.at(i, j, k)->nodeData.centerWeightContribution[0] +=
						(Real)(dx[0][i] * dx[1][j] * weight * dx[2][k]);
			}
		}
	}
}

template<int Degree, bool OutputDensity>
bool Octree<Degree, OutputDensity>::inBounds(Point3D<Real> p) const {
	Real e = boundaryType_ == BoundaryTypeNone ? 0.25 : 0;
	return p[0] >= e && p[0] <= 1 - e && p[1] >= e && p[1] <= 1 - e && p[2] >= e && p[2] <= 1 - e;
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::setTree(std::string const& fileName, int maxDepth, int minDepth,
		int splatDepth, Real samplesPerNode, Real scaleFactor, bool useConfidence,
		bool useNormalWeights, Real constraintWeight, int adaptiveExponent, XForm<Real, 4> xForm) {
	if(splatDepth < 0) splatDepth = 0;
	samplesPerNode_ = samplesPerNode;
	splatDepth_ = splatDepth;
	constrainValues_ = constraintWeight > 0;

	XForm<Real, 3> xFormN = xForm.cut<3>().transpose().inverse();
	if(boundaryType_ == BoundaryTypeNone) {
		++maxDepth;
		minDepth_ = clamp(minDepth + 1, 2, maxDepth);
		if(splatDepth > 0) ++splatDepth;
	} else minDepth_ = clamp(minDepth, 0, maxDepth);

	TreeNeighborKey3 neighborKey(maxDepth);
	PointStream<Real>* pointStream = PointStream<Real>::open(fileName);

	// TODO: PointStream should move to proper c++ iterators
	{
		Point3D<Real> min;
		Point3D<Real> max;
		bool unassigned = true;
		Point3D<Real> p;
		Point3D<Real> n;
		// Read through once to get the center and scale
		while(pointStream->nextPoint(p, n)) {
			p = xForm * p;
			for(int i = 0; i != DIMENSION; ++i) {
				if(unassigned || p[i] < min[i]) min[i] = p[i];
				if(unassigned || p[i] > max[i]) max[i] = p[i];
			}
			unassigned = false;
		}
		scale_ = std::max(max[0] - min[0], std::max(max[1] - min[1], max[2] - min[2]));
		scale_ *= boundaryType_ == BoundaryTypeNone ? 2 * scaleFactor : scaleFactor;
		center_ = (max + min) / 2 - Point3D<Real>::ones() * (scale_ / 2);
	}

	tree_.setFullDepth(minDepth_);
	if(splatDepth > 0) {
		pointStream->reset();
		Point3D<Real> p;
		Point3D<Real> n;
		while(pointStream->nextPoint(p, n)) {
			p = (xForm * p - center_) / scale_;
			n = xFormN * n;
			if(!inBounds(p)) continue;
			Point3D<Real> myCenter(0.5, 0.5, 0.5);
			Real myWidth = 1;
			Real weight = useConfidence ? Length(n) : 1;
			TreeOctNode* temp = &tree_;
			int d = 0;
			while(d < splatDepth) {
				UpdateWeightContribution(temp, p, neighborKey, weight);
				if(!temp->hasChildren()) temp->initChildren();
				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				temp = temp->child(cIndex);
				myWidth /= 2;
				myCenter[0] += (cIndex & 1 ? 1 : -1) * myWidth / 2;
				myCenter[1] += (cIndex & 2 ? 1 : -1) * myWidth / 2;
				myCenter[2] += (cIndex & 4 ? 1 : -1) * myWidth / 2;
				++d;
			}
			UpdateWeightContribution(temp, p, neighborKey, weight);
		}
	}

	double pointWeightSum = 0;
	normals_.clear();
	int cnt = 0;
	pointStream->reset();
	Point3D<Real> p;
	Point3D<Real> n;
	while(pointStream->nextPoint(p, n)) {
		p = (xForm * p - center_) / scale_;
		n = xFormN * (-n);
		if(!inBounds(p)) continue;
		Real normalLength = Length(n);
		if(normalLength <= EPSILON) continue;
		if(!useConfidence) n /= normalLength;

		if(samplesPerNode > 0 && splatDepth) {
			pointWeightSum += SplatOrientedPoint(p, n, neighborKey, splatDepth, samplesPerNode,
					minDepth_, maxDepth);
		} else {
			TreeOctNode* temp = &tree_;
			Point3D<Real> myCenter(0.5, 0.5, 0.5);
			Real myWidth = 1;
			int d = 0;
			if(splatDepth) {
				while(d < splatDepth) {
					int cIndex = TreeOctNode::CornerIndex(myCenter, p);
					temp = temp->child(cIndex);
					myWidth /= 2;
					myCenter[0] += (cIndex & 1 ? 1 : -1) * myWidth / 2;
					myCenter[1] += (cIndex & 2 ? 1 : -1) * myWidth / 2;
					myCenter[2] += (cIndex & 4 ? 1 : -1) * myWidth / 2;
					++d;
				}
				Real pointWeight = 
					GetSampleWeight(temp, p, (TreeConstNeighbors3&)neighborKey.setNeighbors(temp));
				n *= pointWeight;
				pointWeightSum += pointWeight;
			}
			while(d < maxDepth) {
				if(!temp->hasChildren()) temp->initChildren();
				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				temp = temp->child(cIndex);
				myWidth /= 2;
				myCenter[0] += (cIndex & 1 ? 1 : -1) * myWidth / 2;
				myCenter[1] += (cIndex & 2 ? 1 : -1) * myWidth / 2;
				myCenter[2] += (cIndex & 4 ? 1 : -1) * myWidth / 2;
				++d;
			}
			SplatOrientedPoint(temp, p, n, neighborKey);
		}
		if(constrainValues_) {
			Real pointScreeningWeight = useNormalWeights ? normalLength : 1;
			TreeOctNode* temp = &tree_;
			Point3D<Real> myCenter(0.5, 0.5, 0.5);
			Real myWidth = 1;
			int d = 0;
			while(1) {
				int idx = temp->nodeData.pointIndex;
				if(idx == -1) {
					idx = points_.size();
					points_.push_back(PointData(p * pointScreeningWeight, pointScreeningWeight));
					temp->nodeData.pointIndex = idx;
				} else {
					points_[idx].weight += pointScreeningWeight;
					points_[idx].position += p * pointScreeningWeight;
				}

				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				if(!temp->hasChildren()) break;
				temp = temp->child(cIndex);
				myWidth /= 2;
				myCenter[0] += (cIndex & 1 ? 1 : -1) * myWidth / 2;
				myCenter[1] += (cIndex & 2 ? 1 : -1) * myWidth / 2;
				myCenter[2] += (cIndex & 4 ? 1 : -1) * myWidth / 2;
				++d;
			}
		}
		++cnt;
	}

	if(boundaryType_ == BoundaryTypeNone) pointWeightSum *= 4;
	constraintWeight *= pointWeightSum / cnt;

	MemoryUsage();
	delete pointStream;
	if(constrainValues_)
		for(TreeOctNode* node = tree_.nextNode(); node; node = tree_.nextNode(node))
			if(node->nodeData.pointIndex != -1) {
				int idx = node->nodeData.pointIndex;
				points_[idx].position /= points_[idx].weight;
				int nd = boundaryType_ == BoundaryTypeNone ? node->depth() - 1 : node->depth();
				int md = boundaryType_ == BoundaryTypeNone ? maxDepth - 1 : maxDepth;
				int e = nd * adaptiveExponent - md * (adaptiveExponent - 1);
				Real mul = e < 0 ? (Real)1 / (1 << (-e)) : 1 << e;
				points_[idx].weight *= mul * constraintWeight;
			}
#if FORCE_NEUMANN_FIELD
	if(boundaryType_ == BoundaryTypeNeumann)
		for(TreeOctNode* node = tree_.nextNode(); node; node = tree_.nextNode(node)) {
			int d;
			int off[3];
			node->depthAndOffset(d, off);
			int res = 1 << d;
			if(node->nodeData.normalIndex < 0) continue;
			Point3D<Real>& normal = normals_[node->nodeData.normalIndex];
			for(int i = 0; i != 3; ++i)
				if(off[i] == 0 || off[i] == res - 1) normal[i] = 0;
		}
#endif // FORCE_NEUMANN_FIELD
	MemoryUsage();
	return cnt;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::finalize(int subdivideDepth) {
	int maxDepth = tree_.maxDepth();
	TreeNeighborKey3 nKey(maxDepth);
	for(int d = maxDepth; d > 1; --d)
		for(TreeOctNode* node = tree_.nextNode(); node; node = tree_.nextNode(node))
			if(node->depth() == d) {
				TreeNeighbors3& neighbors = nKey.setNeighbors(node->parent()->parent());
				for(int i = 0; i != 3; ++i)
					for(int j = 0; j != 3; ++j)
						for(int k = 0; k != 3; ++k)
							if(neighbors.at(i, j, k) && !neighbors.at(i, j, k)->hasChildren())
								neighbors.at(i, j, k)->initChildren();
			}
	refineBoundary(subdivideDepth);
}

template<int Degree, bool OutputDensity>
double Octree<Degree, OutputDensity>::GetLaplacian(Integrator const& integrator, int d,
		int const off1[], int const off2[], bool childParent) const {
	double vv[] = {
		integrator.dot(d, off1[0], off2[0], false, false, childParent),
		integrator.dot(d, off1[1], off2[1], false, false, childParent),
		integrator.dot(d, off1[2], off2[2], false, false, childParent)
	};
	double dd[] = {
		integrator.dot(d, off1[0], off2[0], true, true, childParent),
		integrator.dot(d, off1[1], off2[1], true, true, childParent),
		integrator.dot(d, off1[2], off2[2], true, true, childParent)
	};
	return dd[0] * vv[1] * vv[2] + vv[0] * dd[1] * vv[2] + vv[0] * vv[1] * dd[2];
}

template<int Degree, bool OutputDensity>
double Octree<Degree, OutputDensity>::GetDivergence1(Integrator const& integrator, int d,
		int const off1[], int const off2[], bool childParent, Point3D<Real> const& normal1) const {
	return Dot(GetDivergence1(integrator, d, off1, off2, childParent), Point3D<double>(normal1));
}

template<int Degree, bool OutputDensity> 
double Octree< Degree , OutputDensity >::GetDivergence2(Integrator const& integrator, int d,
		int const off1[], int const off2[], bool childParent, Point3D<Real> const& normal2) const {
	return Dot(GetDivergence2(integrator, d, off1, off2, childParent), Point3D<double>(normal2));
}

template<int Degree, bool OutputDensity>
Point3D<double> Octree<Degree, OutputDensity>::GetDivergence1(Integrator const& integrator, int d,
		int const off1[], int const off2[], bool childParent) const {
	double vv[] = {
		integrator.dot(d, off1[0], off2[0], false, false, childParent),
		integrator.dot(d, off1[1], off2[1], false, false, childParent),
		integrator.dot(d, off1[2], off2[2], false, false, childParent)
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double vd[] = {
		integrator.dot(d, off1[0], off2[0], false, true, childParent),
		integrator.dot(d, off1[1], off2[1], false, true, childParent),
		integrator.dot(d, off1[2], off2[2], false, true, childParent)
	};
	return Point3D<double>(vd[0] * vv[1] * vv[2], vv[0] * vd[1] * vv[2], vv[0] * vv[1] * vd[2]);
#else
	// Take the dot-product of the divergence of the vector-field with the basis function
	double dv[] = {
		integrator.dot(d, off1[0], off2[0], true, false, childParent),
		integrator.dot(d, off1[1], off2[1], true, false, childParent),
		integrator.dot(d, off1[2], off2[2], true, false, childParent)
	};
	return -Point3D<double>(dv[0] * vv[1] * vv[2], vv[0] * dv[1] * vv[2], vv[0] * vv[1] * dv[2]);
#endif
}

template<int Degree, bool OutputDensity> 
Point3D<double> Octree<Degree, OutputDensity>::GetDivergence2(Integrator const& integrator, int d,
		int const off1[], int const off2[], bool childParent) const {
	double vv[] = {
		integrator.dot(d, off1[0], off2[0], false, false, childParent),
		integrator.dot(d, off1[1], off2[1], false, false, childParent),
		integrator.dot(d, off1[2], off2[2], false, false, childParent)
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double dv[] = {
		integrator.dot(d, off1[0], off2[0], true, false, childParent),
		integrator.dot(d, off1[1], off2[1], true, false, childParent),
		integrator.dot(d, off1[2], off2[2], true, false, childParent)
	};
	return Point3D<double>(dv[0] * vv[1] * vv[2], vv[0] * dv[1] * vv[2], vv[0] * vv[1] * dv[2]);
#else
	// Take the dot-product of the divergence of the vector-field with the basis function
	double vd[] = {
		integrator.dot(d, off1[0], off2[0], false, true, childParent),
		integrator.dot(d, off1[1], off2[1], false, true, childParent),
		integrator.dot(d, off1[2], off2[2], false, true, childParent)
	};
	return -Point3D<double>(vd[0] * vv[1] * vv[2], vv[0] * vd[1] * vv[2], vv[0] * vv[1] * vd[2]);
#endif
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::SetMatrixRowBounds(TreeOctNode const* node, int rDepth,
		int const rOff[3], Range3D& range) const {
	int depth;
	int off[3];
	node->depthAndOffset(depth, off);
	int width = 1 << (depth - rDepth);
	int* starts[] = { &range.xStart, &range.yStart, &range.zStart };
	int* ends[] = { &range.xEnd, &range.yEnd, &range.zEnd };

	for(int i = 0; i != 3; ++i) {
		off[i] -= rOff[i] << (depth - rDepth);
		*starts[i] = off[i] < 0 ? -off[i] : 0;
		*ends[i] = off[i] >= width ? (4 - (off[i] - width)) : 5;
	}
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::GetMatrixRowSize(TreeNeighbors5 const& neighbors5,
		Range3D const& range, bool symmetric) const {
	int count = 0;
	if(symmetric) {
		for(int x = range.xStart; x != 3; ++x) {
			for(int y = range.yStart; y != range.yEnd; ++y) {
				if(x == 2 && y > 2) continue;
				for(int z = range.zStart; z < range.zEnd; ++z) {
					if(x == 2 && y == 2 && z > 2) continue;
					if(neighbors5.at(x, y, z) && neighbors5.at(x, y, z)->nodeData.nodeIndex >= 0)
						++count;
				}
			}
		}
	} else {
		int nodeIndex = neighbors5.at(2, 2, 2)->nodeData.nodeIndex;
		for(int x = range.xStart; x != range.xEnd; ++x)
			for(int y = range.yStart; y != range.yEnd; ++y)
				for(int z = range.zStart; z != range.zEnd; ++z)
					if(neighbors5.at(x, y, z) && neighbors5.at(x, y, z)->nodeData.nodeIndex >= 0 &&
							(!symmetric || neighbors5.at(x, y, z)->nodeData.nodeIndex >= nodeIndex))
						++count;
	}
	return count;
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::SetMatrixRow(TreeNeighbors5 const& neighbors5,
		SparseSymmetricMatrix<MatrixReal>& m, int row, int offset, Integrator const& integrator,
		Stencil<double, 5> const& stencil, Range3D const& range, bool symmetric) const {
	TreeOctNode const* node = neighbors5.at(2, 2, 2);
	int d;
	int off[3];
	node->depthAndOffset(d, off);

	Real pointValues[5][5][5] = {};
	if(constrainValues_) {
		int idx[3] = {
			BinaryNode<double>::CenterIndex(d, off[0]),
			BinaryNode<double>::CenterIndex(d, off[1]),
			BinaryNode<double>::CenterIndex(d, off[2])
		};
		Real diagonal = 0;
		for(int i = 0; i != 3; ++i) {
			for(int j = 0; j != 3; ++j) {
				for(int k = 0; k != 3; ++k) {
					TreeOctNode const* _node = neighbors5.at(i + 1, j + 1, k + 1);
					if(_node && _node->nodeData.pointIndex != -1) {
						Real splineValues[3][3] = {};
						PointData const& pData = points_[_node->nodeData.pointIndex];
						int jdx[] = { i, j, k };
						for(int l = 0; l != 3; ++l) {
							for(int s = 0; s != 3; ++s) {
								int m = idx[l] + jdx[l] - s;
								if(m >= 0 && m < (2 << node->depth()) - 1)
									splineValues[l][s] = fData_.baseBSplines(m, s)(pData.position[l]);
							}
						}
						Real value = splineValues[0][i] * splineValues[1][j] * splineValues[2][k];
						Real weightedValue = value * pData.weight;
						diagonal += value * value * pData.weight;
						for(int s = 0; s != 3; ++s) splineValues[0][s] *= weightedValue;
						for(int ii = 0; ii != 3; ++ii)
							for(int jj = 0; jj != 3; ++jj)
								for(int kk = 0; kk != 3; ++kk)
									pointValues[i + ii][j + jj][k + kk] +=
										splineValues[0][2 - ii] * splineValues[1][2 - jj] * splineValues[2][2 - kk];
					}
				}
			}
		}
		pointValues[2][2][2] = diagonal;
	}

	int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 2 : 2;
	int mx = (1 << d) - mn;
	bool isInterior =
		off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
	int count = 0;
	for(int x = range.xStart; x < (symmetric ? 3 : range.xEnd); ++x) {
		for(int y = range.yStart; y < range.yEnd; ++y) {
			if(x == 2 && y > 2 && symmetric) break;
			for(int z = range.zStart; z < range.zEnd; ++z) {
				if(x == 2 && y == 2 && z > 2 && symmetric) break;
				TreeOctNode const* _node = neighbors5.at(x, y, z);
				if(_node && _node->nodeData.nodeIndex >= 0) {
					Real temp = 0;
					if(isInterior) temp = stencil.at(x, y, z);
					else {
						int _d;
						int _off[3];
						_node->depthAndOffset(_d, _off);
						temp = GetLaplacian(integrator, d, off, _off, false);
					}
					if(constrainValues_) temp += pointValues[x][y][z];
					if(x == 2 && y == 2 && z == 2 && symmetric) temp /= 2;
					if(std::abs(temp) > MATRIX_ENTRY_EPSILON) {
						m.at(row, count).N = _node->nodeData.nodeIndex - offset;
						m.at(row, count).Value = temp;
						++count;
					}
				}
			}
		}
	}
	return count;
}

template<class T, int N, class F>
Stencil<T, N> SetStencil(F const& func) {
	Stencil<T, N> stencil;
	for(int x = 0; x != N; ++x) {
		for(int y = 0; y != N; ++y) {
			for(int z = 0; z != N; ++z) {
				stencil.at(x, y, z) = func(x, y, z);
			}
		}
	}
	return stencil;
}

template<class T, class F>
struct SetStencil1Function {
	SetStencil1Function(F const& f, int i, int j, int k): func(f), i(i), j(j), k(k) { }
	T operator()(int x, int y, int z) const { return func(i, j, k, x, y, z); }
	F const& func;
	int i;
	int j;
	int k;
};

template<class T, int N, class F>
struct SetStencil2Function {
	SetStencil2Function(F const& f): func(f) { }
	Stencil<T, N> operator()(int i, int j, int k) const {
		return SetStencil<T, N, SetStencil1Function<T, F> >(SetStencil1Function<T, F>(func, i, j, k));
	}
	F const& func;
};

template<class T, int N1, int N2, class F>
Stencil<Stencil<T, N2>, N1> SetStencil(F const& func) {
	return SetStencil<Stencil<T, N2>, N1>(SetStencil2Function<T, N2, F>(func));
}

template<class T, class F>
struct SetStencil3Function {
	SetStencil3Function(F const& f, int cx, int cy, int cz): func(f), cx(cx), cy(cy), cz(cz) { }
	T operator()(int i, int j, int k, int x, int y, int z) const
		{ return func(cx, cy, cz, i, j, k, x, y, z); }
	F const& func;
	int cx;
	int cy;
	int cz;
};

template<class T, int N1, int N2, class F>
struct SetStencil4Function {
	SetStencil4Function(F const& f): func(f) { }
	Stencil<Stencil<T, N2>, N1> operator()(int i, int j, int k) const {
		return SetStencil<T, N1, N2, SetStencil3Function<T, F> >(SetStencil3Function<T, F>(func, i, j, k));
	}
	F const& func;
};

template<class T, int N1, int N2, int N3, class F>
Stencil<Stencil<Stencil<T, N3>, N2>, N1> SetStencil(F const& func) {
	return SetStencil<Stencil<Stencil<T, N3>, N2>, N1>(SetStencil4Function<T, N2, N3, F>(func));
}

template<int Degree, bool OutputDensity>
DivergenceStencil Octree<Degree, OutputDensity>::SetDivergenceStencil(int depth,
		Integrator const& integrator, bool scatter) const {
	if(depth < 2) return DivergenceStencil();
	int center = 1 << (depth - 1);
	return SetStencil<Point3D<double>, 5, SetDivergenceStencilFunction>(
		SetDivergenceStencilFunction(*this, depth, integrator, scatter, center));
}

template<int Degree, bool OutputDensity>
DivergenceStencils Octree<Degree, OutputDensity>::SetDivergenceStencils(int depth,
		Integrator const& integrator, bool scatter) const {
	if(depth < 2) return DivergenceStencils();
	int center = 1 << (depth - 1);
	return SetStencil<Point3D<double>, 2, 5, SetDivergenceStencilsFunction>(
		SetDivergenceStencilsFunction(*this, depth, integrator, scatter, center));
}

template<int Degree, bool OutputDensity>
LaplacianStencil Octree<Degree, OutputDensity>::SetLaplacianStencil(int depth,
		Integrator const& integrator) const {
	if(depth < 2) return LaplacianStencil();
	int center = 1 << (depth - 1);
	return SetStencil<double, 5, SetLaplacianStencilFunction>(
		SetLaplacianStencilFunction(*this, depth, integrator, center));
}

template<int Degree, bool OutputDensity>
LaplacianStencils Octree<Degree, OutputDensity>::SetLaplacianStencils(int depth,
		Integrator const& integrator) const {
	if(depth < 2) return LaplacianStencils();
	int center = 1 << (depth - 1);
	return SetStencil<double, 2, 5, SetLaplacianStencilsFunction>(
		SetLaplacianStencilsFunction(*this, depth, integrator, center));
}

template<int Degree, bool OutputDensity>
CenterEvaluationStencil Octree<Degree, OutputDensity>::SetCenterEvaluationStencil(
		CenterEvaluator1 const& evaluator, int depth) const {
	if(depth < 2) return CenterEvaluationStencil();
	int center = 1 << (depth - 1);
	return SetStencil<double, 3, SetCenterEvaluationStencilFunction>(
		SetCenterEvaluationStencilFunction(depth, center, evaluator));
}

template<int Degree, bool OutputDensity>
CenterEvaluationStencils Octree<Degree, OutputDensity>::SetCenterEvaluationStencils(
		CenterEvaluator1 const& evaluator, int depth) const {
	if(depth < 3) return CenterEvaluationStencils();
	int center = 1 << (depth - 1);
	return SetStencil<double, 2, 3, SetCenterEvaluationStencilsFunction>(
		SetCenterEvaluationStencilsFunction(depth, center, evaluator));
}

template<int Degree, bool OutputDensity>
CornerEvaluationStencil Octree<Degree, OutputDensity>::SetCornerEvaluationStencil(
		CornerEvaluator2 const& evaluator, int depth) const {
	if(depth < 2) return CornerEvaluationStencil();
	int center = 1 << (depth - 1);
	return SetStencil<double, 2, 3, SetCornerEvaluationStencilFunction>(
		SetCornerEvaluationStencilFunction(depth, center, evaluator));
}

template<int Degree, bool OutputDensity>
CornerEvaluationStencils Octree<Degree, OutputDensity>::SetCornerEvaluationStencils(
		CornerEvaluator2 const& evaluator, int depth) const {
	if(depth < 3) return CornerEvaluationStencils();
	int center = 1 << (depth - 1);
	return SetStencil<double, 2, 2, 3, SetCornerEvaluationStencilsFunction>(
		SetCornerEvaluationStencilsFunction(depth, center, evaluator));
}

template<int Degree, bool OutputDensity>
CornerNormalEvaluationStencil Octree<Degree, OutputDensity>::SetCornerNormalEvaluationStencil(
		CornerEvaluator2 const& evaluator, int depth) const {
	if(depth < 2) return CornerNormalEvaluationStencil();
	int center = 1 << (depth - 1);
	return SetStencil<Point3D<double>, 2, 5, SetCornerNormalEvaluationStencilFunction>(
		SetCornerNormalEvaluationStencilFunction(depth, center, evaluator));
}

template<int Degree, bool OutputDensity>
CornerNormalEvaluationStencils Octree<Degree, OutputDensity>::SetCornerNormalEvaluationStencils(
		CornerEvaluator2 const& evaluator, int depth) const {
	if(depth < 3) return CornerNormalEvaluationStencils();
	int center = 1 << (depth - 1);
	return SetStencil<Point3D<double>, 2, 2, 5, SetCornerNormalEvaluationStencilsFunction>(
			SetCornerNormalEvaluationStencilsFunction(depth, center, evaluator));
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::UpdateCoarserSupportBounds(TreeOctNode const* node, Range3D& range) {
	if(!node->parent()) return;
	int x;
	int y;
	int z;
	Cube::FactorCornerIndex(node->parent()->childIndex(node), x, y, z);
	if(x == 0) range.xEnd = 4;
	else range.xStart = 1;
	if(y == 0) range.yEnd = 4;
	else range.yStart = 1;
	if(z == 0) range.zEnd = 4;
	else range.zStart = 1;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::UpdateConstraintsFromCoarser(TreeNeighbors5 const& neighbors5,
		TreeNeighbors5 const& pNeighbors5, TreeOctNode* node, Real const* metSolution,
		Integrator const& integrator, Stencil<double, 5> const& lapStencil) const {
	int d;
	int off[3];
	node->depthAndOffset(d, off);
	int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 4 : 4;
	int mx = (1 << d) - mn;
	bool isInterior =
		off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
	if(d <= minDepth_) return;
	// Offset the constraints using the solution from lower resolutions.
	Range3D range = Range3D::FullRange();
	UpdateCoarserSupportBounds(node, range);

	for(int x = range.xStart; x != range.xEnd; ++x) {
		for(int y = range.yStart; y != range.yEnd; ++y) {
			for(int z = range.zStart; z != range.zEnd; ++z) {
				TreeOctNode const* _node = pNeighbors5.at(x, y, z);
				if(_node && _node->nodeData.nodeIndex >= 0) {
					Real _solution = metSolution[_node->nodeData.nodeIndex];
					if(isInterior) node->nodeData.constraint -= (Real)(lapStencil.at(x, y, z) * _solution);
					else {
						int _d;
						int _off[3];
						_node->depthAndOffset(_d, _off);
						node->nodeData.constraint -=
							(Real)(GetLaplacian(integrator, d, off, _off, true) * _solution);
					}
				}
			}
		}
	}
	if(constrainValues_) {
		double constraint = 0;
		off[0] = BinaryNode<double>::CenterIndex(d, off[0]);
		off[1] = BinaryNode<double>::CenterIndex(d, off[1]);
		off[2] = BinaryNode<double>::CenterIndex(d, off[2]);
		for(int x = 1; x != 4; ++x) {
			for(int y = 1; y != 4; ++y) {
				for(int z = 1; z != 4; ++z) {
					TreeOctNode const* _node = neighbors5.at(x, y, z);
					if(_node && _node->nodeData.pointIndex != -1) {
						PointData const& pData = points_[_node->nodeData.pointIndex];
						Real pointValue = pData.coarserValue;
						Point3D<Real> p = pData.position;
						constraint += fData_.baseBSplines(off[0], x - 1)(p[0]) *
							fData_.baseBSplines(off[1], y - 1)(p[1]) *
							fData_.baseBSplines(off[2], z - 1)(p[2]) *
							pointValue;
					}
				}
			}
		}
		node->nodeData.constraint -= (Real)constraint;
	}
}

template<bool OutputDensity, class TreeOctNode, class F>
void UpSampleGeneric(int depth, SortedTreeNodes<OutputDensity> const& sNodes, BoundaryType boundaryType,
		int threads, F const& func) {
	double cornerValue = boundaryType == BoundaryTypeDirichlet ? 0.5 :
		boundaryType == BoundaryTypeNeumann ? 1 : 0.75;
	// For every node at the current depth
	typename TreeOctNode::NeighborKey3 neighborKey(depth);
#pragma omp parallel for num_threads(threads) firstprivate(neighborKey)
	for(int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; ++i) {
		int d;
		int off[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		UpSampleData usData[3];
		for(int dd = 0; dd != 3; ++dd)
			usData[dd] = off[dd] == 0 ? UpSampleData(1, cornerValue, 0) :
				off[dd] + 1 == (1 << depth) ? UpSampleData(0, 0, cornerValue) :
				off[dd] % 2 ? UpSampleData(1, 0.75, 0.25) :
				UpSampleData(0, 0.25, 0.75);
		typename TreeOctNode::Neighbors3& neighbors =
			neighborKey.getNeighbors3(sNodes.treeNodes[i]->parent());
		for(int ii = 0; ii != 2; ++ii) {
			for(int jj = 0; jj != 2; ++jj) {
				for(int kk = 0; kk != 2; ++kk) {
					TreeOctNode const* node =
						neighbors.at(ii + usData[0].start, jj + usData[1].start, kk + usData[2].start);
					if(node && node->nodeData.nodeIndex != -1) {
						int idx[] = { ii, jj, kk };
						func(i, node, usData, idx);
					}
				}
			}
		}
	}
}

template<int Degree, bool OutputDensity>
Vector<Real> Octree<Degree, OutputDensity>::UpSampleCoarserSolution(int depth,
		SortedTreeNodes<OutputDensity> const& sNodes) const {
	size_t start = sNodes.nodeCount[depth];
	size_t end = sNodes.nodeCount[depth + 1];
	Vector<Real> Solution(end - start);
	if((boundaryType_ != BoundaryTypeNone && depth == 0) ||
			(boundaryType_ == BoundaryTypeNone && depth <= 2)) return Solution;
	UpSampleGeneric<OutputDensity, TreeOctNode>(depth, sNodes, boundaryType_, threads_,
		UpSampleCoarserSolutionFunction(Solution, start));
	// Clear the coarser solution
#pragma omp parallel for num_threads(threads_)
	for(int i = sNodes.nodeCount[depth - 1]; i < (int)sNodes.nodeCount[depth]; ++i)
		sNodes.treeNodes[i]->nodeData.solution = 0;
	return Solution;
}

template<int Degree, bool OutputDensity>
template<class C>
void Octree<Degree, OutputDensity>::DownSample(int depth, SortedTreeNodes<OutputDensity> const& sNodes,
		C* constraints) const {
	if(depth == 0) return;
	UpSampleGeneric<OutputDensity, TreeOctNode>(depth, sNodes, boundaryType_, threads_,
		DownSampleFunction<C>(constraints));
}

template<int Degree, bool OutputDensity>
template<class C>
void Octree<Degree, OutputDensity>::UpSample(int depth, SortedTreeNodes<OutputDensity> const& sNodes,
		C* coefficients) const {
	if((boundaryType_ != BoundaryTypeNone && depth == 0) ||
			(boundaryType_ == BoundaryTypeNone && depth <= 2)) return;
	UpSampleGeneric<OutputDensity, TreeOctNode>(depth, sNodes, boundaryType_, threads_,
		UpSample1Function<C>(coefficients));
}

template<int Degree, bool OutputDensity>
template<class C>
void Octree<Degree, OutputDensity>::UpSample(int depth, SortedTreeNodes<OutputDensity> const& sNodes,
		C const* coarseCoefficients, C* fineCoefficients) const {
	if(depth <= minDepth_) return;
	UpSampleGeneric<OutputDensity, TreeOctNode>(depth, sNodes, boundaryType_, threads_,
		UpSample2Function<C>(depth, fineCoefficients, coarseCoefficients, sNodes));
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::SetCoarserPointValues(int depth,
		SortedTreeNodes<OutputDensity> const& sNodes, Real* metSolution) {
	TreeNeighborKey3 neighborKey(depth);
#pragma omp parallel for num_threads(threads_) firstprivate(neighborKey)
	for(int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; ++i) {
		TreeOctNode* node = sNodes.treeNodes[i];
		if(node->nodeData.pointIndex != -1) {
			neighborKey.getNeighbors3(node);
			points_[node->nodeData.pointIndex].coarserValue =
				WeightedCoarserFunctionValue(neighborKey, node, metSolution);
		}
	}
}

template<int Degree, bool OutputDensity>
Real Octree<Degree, OutputDensity>::WeightedCoarserFunctionValue(TreeNeighborKey3 const& neighborKey,
		TreeOctNode const* pointNode, Real* metSolution) const {
	int depth = pointNode->depth();
	if(boundaryType_ == BoundaryTypeDirichlet && depth == 0 && pointNode->nodeData.pointIndex != -1)
		return (Real)(-0.5) * points_[pointNode->nodeData.pointIndex].weight;
	if((boundaryType_ != BoundaryTypeNone && depth == 0) ||
			(boundaryType_ == BoundaryTypeNone && depth <= 2) || pointNode->nodeData.pointIndex == -1)
		return 0;

	Real weight = points_[pointNode->nodeData.pointIndex].weight;
	Point3D<Real> p = points_[pointNode->nodeData.pointIndex].position;
	double pointValue = 0;

	// Iterate over all basis functions that overlap the point at the coarser resolutions
	int d;
	int _idx[3];
	TreeNeighbors3 const& neighbors = neighborKey.neighbors(depth - 1);
	neighbors.at(1, 1, 1)->depthAndOffset(d, _idx);
	_idx[0] = BinaryNode<double>::CenterIndex(d, _idx[0] - 1);
	_idx[1] = BinaryNode<double>::CenterIndex(d, _idx[1] - 1);
	_idx[2] = BinaryNode<double>::CenterIndex(d, _idx[2] - 1);

	for(int j = 0; j != 3; ++j) {
		if(!(_idx[0] + j >= 0 && _idx[0] + j < (1 << depth) - 1)) continue;
		double xValue = fData_.baseBSplines(_idx[0] + j, 2 - j)(p[0]);
		for(int k = 0; k != 3; ++k) {
			if(!(_idx[1] + k >= 0 && _idx[1] + k < (1 << depth) - 1)) continue;
			double xyValue = xValue * fData_.baseBSplines(_idx[1] + k, 2 - k)(p[1]);
			double _pointValue = 0;
			for(int l = 0; l != 3; ++l) {
				TreeOctNode const* basisNode = neighbors.at(j, k, l);
				if(basisNode && basisNode->nodeData.nodeIndex >= 0 &&
						_idx[2] + l >= 0 && _idx[2] + l < (1 << depth) - 1)
					_pointValue += fData_.baseBSplines(_idx[2] + l, 2 - l)(p[2]) *
						(double)metSolution[basisNode->nodeData.nodeIndex];
			}
			pointValue += _pointValue * xyValue;
		}
	}
	if(boundaryType_ == BoundaryTypeDirichlet) pointValue -= (Real)0.5;
	return pointValue * weight;
}

template<int Degree, bool OutputDensity>
template<class F1, class F2, class F3>
SparseSymmetricMatrix<Real> Octree<Degree, OutputDensity>::GetFixedDepthLaplacianGeneric(int depth,
		Integrator const& integrator, SortedTreeNodes<OutputDensity> const& sNodes,
		Real const* metSolution, size_t range, F1 const& getNode, F2 const& getRowSize, F3 const& setRow) {
	SparseSymmetricMatrix<Real> matrix;
	matrix.Resize(range);
	Stencil<double, 5> stencil = SetLaplacianStencil(depth, integrator);
	LaplacianStencils stencils = SetLaplacianStencils(depth, integrator);
	TreeNeighborKey3 neighborKey3(depth);
#pragma omp parallel for num_threads(threads_) firstprivate(neighborKey3)
	for(int i = 0; i < (int)range; ++i) {
		TreeOctNode* node = getNode(i);

		// Get the matrix row size
		bool insetSupported = boundaryType_ != BoundaryTypeNone || IsInsetSupported(node);
		TreeNeighbors5 neighbors5;
		if(insetSupported) neighbors5 = neighborKey3.getNeighbors5(node);
		int count = insetSupported ? getRowSize(neighbors5, true) : 1;

		// Allocate memory for the row
#pragma omp critical(matrix_set_row_size)
		matrix.SetRowSize(i, count);

		// Set the row entries
		if(insetSupported) matrix.rowSize(i) =
				setRow(neighbors5, matrix, i, sNodes.nodeCount[depth], integrator, stencil, true);
		else {
			matrix.at(i, 0) = MatrixEntry<Real>(i, 1);
			matrix.rowSize(i) = 1;
		}

		// Offset the constraints using the solution from lower resolutions.
		int x = 0;
		int y = 0;
		int z = 0;
		if(node->parent()) {
			int c = node->parent()->childIndex(node);
			Cube::FactorCornerIndex(c, x, y, z);
		}
		if(insetSupported) {
			TreeNeighbors5 pNeighbors5 = neighborKey3.getNeighbors5(node->parent());
			UpdateConstraintsFromCoarser(neighbors5, pNeighbors5, node, metSolution, integrator,
					stencils.at(x, y, z));
		}
	}
	return matrix;
}

template<int Degree, bool OutputDensity>
SparseSymmetricMatrix<Real> Octree<Degree, OutputDensity>::GetFixedDepthLaplacian( int depth,
		Integrator const& integrator, SortedTreeNodes<OutputDensity> const& sNodes,
		Real const* metSolution) {
	size_t start = sNodes.nodeCount[depth];
	size_t end = sNodes.nodeCount[depth + 1];
	size_t range = end - start;
	return GetFixedDepthLaplacianGeneric(depth, integrator, sNodes, metSolution, range,
			GetFixedDepthLaplacianGetNodeFunction(sNodes, start),
			GetFixedDepthLaplacianGetRowSizeFunction(*this),
			GetFixedDepthLaplacianSetRowFunction(*this));
}

template<int Degree, bool OutputDensity>
SparseSymmetricMatrix<Real> Octree<Degree, OutputDensity>::GetRestrictedFixedDepthLaplacian(int depth,
		Integrator const& integrator, int const* entries, int entryCount, TreeOctNode const* rNode, Real,
		SortedTreeNodes<OutputDensity> const& sNodes, Real const* metSolution) {
	for(int i = 0; i != (int)entryCount; ++i) sNodes.treeNodes[entries[i]]->nodeData.nodeIndex = i;
	int rDepth;
	int rOff[3];
	rNode->depthAndOffset(rDepth, rOff);
	Range3D range;
	SparseSymmetricMatrix<Real> matrix = GetFixedDepthLaplacianGeneric(depth, integrator, sNodes,
			metSolution, entryCount,
			GetRestrictedFixedDepthLaplacianGetNodeFunction(*this, sNodes, depth, entries,
				rDepth, rOff, range),
			GetRestrictedFixedDepthLaplacianGetRowSizeFunction(*this, range),
			GetRestrictedFixedDepthLaplacianSetRowFunction(*this, range));
	for(int i = 0; i != (int)entryCount; ++i) sNodes.treeNodes[entries[i]]->nodeData.nodeIndex = entries[i];
	return matrix;
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::LaplacianMatrixIteration(int subdivideDepth, bool showResidual,
		int minIters, double accuracy, int maxSolveDepth, int fixedIters) {
	int iter = 0;
	Integrator integrator;
	fData_.setIntegrator(integrator, boundaryType_ == BoundaryTypeNone);
	if(boundaryType_ == BoundaryTypeNone) {
		++subdivideDepth;
		++maxSolveDepth;
	}

	sNodes_.treeNodes[0]->nodeData.solution = 0;

	std::vector<Real> metSolution(sNodes_.nodeCount[sNodes_.maxDepth], 0);
	for(int d = (boundaryType_ == BoundaryTypeNone ? 2 : 0); d != sNodes_.maxDepth; ++d) {
		DumpOutput::instance()("#Depth[%d/%d]: %d\n", boundaryType_ == BoundaryTypeNone ? d - 1 : d,
				boundaryType_ == BoundaryTypeNone ? sNodes_.maxDepth - 2 : sNodes_.maxDepth - 1,
				sNodes_.nodeCount[d + 1] - sNodes_.nodeCount[d]);
		if(subdivideDepth > 0)
			iter += SolveFixedDepthMatrix(d, integrator, sNodes_, &metSolution[0], subdivideDepth,
					showResidual, minIters, accuracy, d > maxSolveDepth, fixedIters);
		else
			iter += SolveFixedDepthMatrix(d, integrator, sNodes_, &metSolution[0],
					showResidual, minIters, accuracy, d > maxSolveDepth, fixedIters);
	}
	return iter;
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::SolveFixedDepthMatrix(int depth, Integrator const& integrator,
		SortedTreeNodes<OutputDensity> const& sNodes, Real* metSolution, bool showResidual, int minIters,
		double accuracy, bool noSolve, int fixedIters) {
	Vector<Real> X(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	if(depth <= minDepth_) X = UpSampleCoarserSolution(depth, sNodes);
	else {
		// Up-sample the cumulative solution into the previous depth
		UpSample(depth - 1, sNodes, metSolution);
		// Add in the solution from that depth
		if(depth)
#pragma omp parallel for num_threads(threads_)
			for(int i = sNodes_.nodeCount[depth - 1]; i < sNodes_.nodeCount[depth]; ++i)
				metSolution[i] += sNodes_.treeNodes[i]->nodeData.solution;
	}
	double evaluateTime = 0;
	if(constrainValues_) {
		evaluateTime = Time();
		SetCoarserPointValues(depth, sNodes, metSolution);
		evaluateTime = Time() - evaluateTime;
	}

	double systemTime = Time();
	// Get the system matrix
	SparseSymmetricMatrix<Real> M = GetFixedDepthLaplacian(depth, integrator, sNodes, metSolution);
	// Set the constraint vector
	Vector<Real> B(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	for(int i = sNodes.nodeCount[depth]; i != sNodes.nodeCount[depth + 1]; ++i)
		B[i - sNodes.nodeCount[depth]] =
			boundaryType_ != BoundaryTypeNone || IsInsetSupported(sNodes.treeNodes[i]) ?
			sNodes.treeNodes[i]->nodeData.constraint : 0;
	systemTime = Time() - systemTime;

	double solveTime = Time();
	// Solve the linear system
	Real _accuracy = (Real)(accuracy / 100000) * M.Rows();
	int res = 1 << depth;
	if(boundaryType_ == BoundaryTypeNone && depth > 3) res -= 1 << (depth - 2);
	int iter = 0;
	if(!noSolve) {
		int iters = fixedIters >= 0 ? fixedIters :
			std::max((int)std::pow(M.Rows(), ITERATION_POWER), minIters);
		Real accuracy = fixedIters >= 0 ? 1e-10 : _accuracy;
		iter += SparseSymmetricMatrix<Real>::Solve(M, B, iters, X, accuracy, false, threads_,
			M.Rows() == res * res * res && !constrainValues_ && boundaryType_ != BoundaryTypeDirichlet);
	}
	solveTime = Time() - solveTime;

	if(showResidual) {
		double bNorm = B.Norm(2);
		double rNorm = (B - M * X).Norm(2);
		DumpOutput::instance()("#\tResidual: (%d %g) %g -> %g (%f) [%d]\n", M.Entries(),
				std::sqrt(M.Norm(2)), bNorm, rNorm, rNorm / bNorm, iter);
	}

	// Copy the solution back into the tree (over-writing the constraints)
	for(int i = sNodes.nodeCount[depth]; i != sNodes.nodeCount[depth+1]; ++i)
		sNodes.treeNodes[i]->nodeData.solution = X[i - sNodes.nodeCount[depth]];

	DumpOutput::instance()("#\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n",
			evaluateTime, systemTime, solveTime, (float)MemoryUsage());
	return iter;
}

template<class TreeOctNode, class ExtraCondition, class DoWork>
void getAdjacencyCount(TreeOctNode* node, typename TreeOctNode::NeighborKey3& neighborKey3, int depth,
		int fDataDepth, int width, ExtraCondition const& extraCondition, DoWork const& doWork) {
	// Count the number of nodes at depth "depth" that lie under sNodes.treeNodes[i]
	for(TreeOctNode* temp = node->nextNode(); temp;) {
		if(temp->depth() == depth && extraCondition(temp)) {
			doWork(temp, nullptr);
			temp = node->nextBranch(temp);
		} else temp = node->nextNode(temp);
	}
	// [ERROR] Wow!!!! This is as stupid as stupid gets. All pairs? Really?
#pragma message("[WARNING] Assuming that the 2-ring contains all the children() of interest...")
	typename TreeOctNode::Neighbors5 neighbors5 = neighborKey3.getNeighbors5(node);
	for(int x = 0; x != 5; ++x)
		for(int y = 0; y != 5; ++y)
			for(int z = 0; z != 5; ++z)
				if(neighbors5.at(x, y, z) && !(x == 2 && y == 2 && z == 2))
					TreeOctNode::ProcessFixedDepthNodeAdjacentNodes(fDataDepth, node,
							1, neighbors5.at(x, y, z), 2 * width - 1, depth, doWork);
}

template<class TreeOctNode>
bool SolveFixedDepthMatrix1Function(TreeOctNode const*) { return true; }

template<class TreeOctNode>
class SolveFixedDepthMatrix2Function {
public:
	SolveFixedDepthMatrix2Function(int& adjacencyCount): adjacencyCount(adjacencyCount) {}
	void operator()(TreeOctNode const*, TreeOctNode const*) const {
		++adjacencyCount;
	}
private:
	int& adjacencyCount;
};

template<class TreeOctNode>
bool SolveFixedDepthMatrix3Function(TreeOctNode const* temp) { return temp->nodeData.nodeIndex != -1; }

template<class TreeOctNode>
class SolveFixedDepthMatrix4Function {
public:
	SolveFixedDepthMatrix4Function(int& adjacencyCount2, int* adjacencies):
		adjacencyCount2(adjacencyCount2), adjacencies(adjacencies) { }
	void operator()(TreeOctNode const* node1, TreeOctNode const*) const {
		adjacencies[adjacencyCount2++] = node1->nodeData.nodeIndex;
	}
private:
	int& adjacencyCount2;
	int* adjacencies;
};

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::SolveFixedDepthMatrix(int depth, Integrator const& integrator,
		SortedTreeNodes<OutputDensity> const& sNodes, Real* metSolution, int startingDepth,
		bool showResidual, int minIters, double accuracy, bool noSolve, int fixedIters) {
	if(startingDepth >= depth)
		return SolveFixedDepthMatrix(depth, integrator, sNodes, metSolution, showResidual, minIters,
				accuracy, noSolve, fixedIters);

	if(depth > minDepth_) {
		// Up-sample the cumulative solution into the previous depth
		UpSample(depth - 1, sNodes, metSolution);
		// Add in the solution from that depth
		if(depth)
#pragma omp parallel for num_threads(threads_)
			for(int i = sNodes_.nodeCount[depth - 1]; i < sNodes_.nodeCount[depth]; ++i)
				metSolution[i] += sNodes_.treeNodes[i]->nodeData.solution;
	}

	double evaluateTime = 0;
	if(constrainValues_) {
		evaluateTime = Time();
		SetCoarserPointValues(depth, sNodes, metSolution);
		evaluateTime = Time() - evaluateTime;
	}

	Vector<Real> B(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	// Back-up the constraints
	for(int i = sNodes.nodeCount[depth]; i != sNodes.nodeCount[depth + 1]; ++i) {
		B[i - sNodes.nodeCount[depth]] =
			boundaryType_ != BoundaryTypeNone || IsInsetSupported(sNodes.treeNodes[i]) ?
			sNodes.treeNodes[i]->nodeData.constraint : 0;
		sNodes.treeNodes[i]->nodeData.constraint = 0;
	}

	int d = depth - startingDepth;
	if(boundaryType_ == BoundaryTypeNone) ++d;
	std::vector<int> subDimension;
	int maxDimension = 0;
	TreeNeighborKey3 neighborKey3(fData_.depth());
	for(int i = sNodes.nodeCount[d]; i != sNodes.nodeCount[d + 1]; ++i) {
		int adjacencyCount = 0;
		getAdjacencyCount<TreeOctNode>(sNodes.treeNodes[i], neighborKey3, depth, fData_.depth(), width_,
				SolveFixedDepthMatrix1Function<TreeOctNode>,
				SolveFixedDepthMatrix2Function<TreeOctNode>(adjacencyCount));
		subDimension.push_back(adjacencyCount);
		maxDimension = std::max(maxDimension, adjacencyCount);
	}

	Real myRadius = lrint(2 * radius_ - (Real)0.5 - ROUND_EPS) + ROUND_EPS;
	int* adjacencies = new int[maxDimension];
	int tIter = 0;
	double systemTime = 0;
	double solveTime = 0;
	// Iterate through the coarse-level nodes
	for(int i = sNodes.nodeCount[d]; i != sNodes.nodeCount[d + 1]; ++i) {
		// Count the number of nodes at depth "depth" that lie under sNodes.treeNodes[i]
		if(!subDimension[i - sNodes.nodeCount[d]]) continue;
		int iter = 0;
		double time = Time();

		// Set the indices for the nodes under, or near, sNodes.treeNodes[i].
		int adjacencyCount2 = 0;
		getAdjacencyCount<TreeOctNode>(sNodes.treeNodes[i], neighborKey3, depth, fData_.depth(), width_,
				SolveFixedDepthMatrix3Function<TreeOctNode>,
				SolveFixedDepthMatrix4Function<TreeOctNode>(adjacencyCount2, adjacencies));
		// Get the associated constraint vector
		Vector<Real> _B(adjacencyCount2);
		Vector<Real> _X(adjacencyCount2);
#pragma omp parallel for num_threads(threads_) schedule(static)
		for(int j = 0; j < adjacencyCount2; ++j) {
			_B[j] = B[adjacencies[j] - sNodes.nodeCount[depth]];
			_X[j] = sNodes.treeNodes[adjacencies[j]]->nodeData.solution;
		}

		// Get the associated matrix
		SparseSymmetricMatrix<Real> _M = GetRestrictedFixedDepthLaplacian(depth, integrator,
				adjacencies, adjacencyCount2, sNodes.treeNodes[i], myRadius, sNodes, metSolution);
#pragma omp parallel for num_threads(threads_) schedule(static)
		for(int j = 0; j < adjacencyCount2; ++j) {
			_B[j] += sNodes.treeNodes[adjacencies[j]]->nodeData.constraint;
			sNodes.treeNodes[adjacencies[j]]->nodeData.constraint = 0;
		}
		systemTime += Time() - time;

		// Solve the matrix
		// Since we don't have the full matrix, the system shouldn't be singular, so we shouldn't have 
		// to correct it
		time = Time();
		Real _accuracy = (Real)(accuracy / 100000) * _M.Rows();
		if(!noSolve) {
			int iters = fixedIters >= 0 ? fixedIters :
				std::max((int)std::pow(_M.Rows(), ITERATION_POWER), minIters);
			Real accuracy = fixedIters >= 0 ? 1e-10 : _accuracy;
			iter += SparseSymmetricMatrix<Real>::Solve(_M, _B, iters, _X, accuracy, false, threads_, false);
		}
		solveTime += Time() - time;

		if(showResidual) {
			double bNorm = _B.Norm(2);
			double rNorm = (_B - _M * _X).Norm(2);
			DumpOutput::instance()("#\t\tResidual: (%d %g) %g -> %g (%f) [%d]\n", _M.Entries(),
					_M.Norm(2), bNorm, rNorm, rNorm / bNorm, iter);
		}

		// Update the solution for all nodes in the sub-tree
#pragma omp parallel for num_threads(threads_)
		for(int j = 0; j < adjacencyCount2; ++j) {
			TreeOctNode* temp = sNodes.treeNodes[adjacencies[j]];
			while(temp->depth() > sNodes.treeNodes[i]->depth()) temp = temp->parent();
			if(temp->nodeData.nodeIndex >= sNodes.treeNodes[i]->nodeData.nodeIndex)
				sNodes.treeNodes[adjacencies[j]]->nodeData.solution = _X[j];
		}
		MemoryUsage();
		tIter += iter;
	}
	delete[] adjacencies;
	MemoryUsage();
	DumpOutput::instance()("#\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n",
			evaluateTime, systemTime, solveTime, (float)maxMemoryUsage());
	return tIter;
}

template<int Degree, bool OutputDensity>
bool Octree<Degree, OutputDensity>::HasNormals(TreeOctNode* node, Real epsilon) const {
	if(node->nodeData.normalIndex >= 0 && (normals_[node->nodeData.normalIndex][0] != 0 ||
				normals_[node->nodeData.normalIndex][1] != 0 ||
				normals_[node->nodeData.normalIndex][2] != 0)) return true;
	if(!node->hasChildren()) return false;
	for(unsigned i = 0; i != Cube::CORNERS; ++i)
		if(HasNormals(node->child(i), epsilon)) return true;
	return false;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::ClipTree() {
	int maxDepth = tree_.maxDepth();
	for(TreeOctNode* temp = tree_.nextNode(); temp; temp = tree_.nextNode(temp)) {
		if(!temp->hasChildren() || temp->depth() < minDepth_) continue;
		bool hasNormals = false;
		for(unsigned i = 0; i != Cube::CORNERS; ++i) {
			if(HasNormals(temp->child(i), EPSILON / (1 << maxDepth))) {
				hasNormals = true;
				break;
			}
		}
		if(!hasNormals) temp->nullChildren();
	}
	MemoryUsage();
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::SetLaplacianConstraints() {
	// To set the Laplacian constraints, we iterate over the
	// splatted normals and compute the dot-product of the
	// divergence of the normal field with all the basis functions.
	// Within the same depth: set directly as a gather
	// Coarser depths 
	Integrator integrator;
	fData_.setIntegrator(integrator, boundaryType_ == BoundaryTypeNone);
	int maxDepth = sNodes_.maxDepth - 1;
	std::vector<Real> constraints(sNodes_.nodeCount[maxDepth]);

	// Clear the constraints
#pragma omp parallel for num_threads(threads_)
	for(int i = 0; i < sNodes_.nodeCount[maxDepth + 1]; ++i)
		sNodes_.treeNodes[i]->nodeData.constraint = 0;

	for(int d = maxDepth; d >= (boundaryType_ == BoundaryTypeNone ? 2 : 0); --d) {
		DivergenceStencil stencil = SetDivergenceStencil(d, integrator, false);
		DivergenceStencils stencils = SetDivergenceStencils(d, integrator, true);
		TreeNeighborKey3 neighborKey3(fData_.depth());
#pragma omp parallel for num_threads(threads_) firstprivate(neighborKey3)
		for(int i = sNodes_.nodeCount[d]; i < sNodes_.nodeCount[d + 1]; ++i) {
			TreeOctNode* node = sNodes_.treeNodes[i];
			Range3D range = Range3D::FullRange();
			TreeNeighbors5 neighbors5 = neighborKey3.getNeighbors5(node);

			int off[3];
			node->depthAndOffset(d, off);
			int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 2 : 2;
			int mx = (1 << d) - mn;
			bool isInterior =
				off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
			mn += 2;
			mx -= 2;
			bool isInterior2 =
				off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
			int cx = 0;
			int cy = 0;
			int cz = 0;
			if(d)
				Cube::FactorCornerIndex(node->parent()->childIndex(node), cx, cy, cz);
			DivergenceStencil& _stencil = stencils.at(cx, cy, cz);

			// Set constraints from current depth
			for(int x = range.xStart; x != range.xEnd; ++x) {
				for(int y = range.yStart; y != range.yEnd; ++y) {
					for(int z = range.zStart; z != range.zEnd; ++z) {
						TreeOctNode const* _node = neighbors5.at(x, y, z);
						if(_node && _node->nodeData.normalIndex >= 0) {
							Point3D<Real> const& _normal = normals_[_node->nodeData.normalIndex];
							int _d;
							int _off[3];
							_node->depthAndOffset(_d, _off);
							node->nodeData.constraint += isInterior ?
								(Real)Dot(stencil.at(x, y, z), Point3D<double>(_normal)) :
								(Real)GetDivergence2(integrator, d, off, _off, false, _normal);
						}
					}
				}
			}
			UpdateCoarserSupportBounds(neighbors5.at(2, 2, 2), range);
			if(node->nodeData.nodeIndex < 0 || node->nodeData.normalIndex < 0) continue;
			Point3D<Real> const& normal = normals_[node->nodeData.normalIndex];
			if(normal == Point3D<Real>()) continue;

			// Set the constraints for the parents
			if(d) {
				neighbors5 = neighborKey3.getNeighbors5(node->parent());

				for(int x = range.xStart; x != range.xEnd; ++x) {
					for(int y = range.yStart; y != range.yEnd; ++y) {
						for(int z = range.zStart; z != range.zEnd; ++z) {
							TreeOctNode* _node = neighbors5.at(x, y, z);
							if(_node && _node->nodeData.nodeIndex != -1) {
								int _d;
								int _off[3];
								_node->depthAndOffset(_d, _off);
								Real c = isInterior2 ?
									Dot(_stencil.at(x, y, z), Point3D<double>(normal)) :
									GetDivergence1(integrator, d, off, _off, true, normal);
#pragma omp atomic
								constraints[_node->nodeData.nodeIndex] += c;
							}
						}
					}
				}
			}
		}
	}
	std::vector<Point3D<Real> > coefficients(sNodes_.nodeCount[maxDepth], Point3D<Real>());
	for(int d = maxDepth - 1; d >= 0; --d) {
#pragma omp parallel for num_threads(threads_)
		for(int i = sNodes_.nodeCount[d]; i < sNodes_.nodeCount[d + 1]; ++i) {
			TreeOctNode* node = sNodes_.treeNodes[i];
			if(node->nodeData.nodeIndex < 0 || node->nodeData.normalIndex < 0) continue;
			coefficients[i] += normals_[node->nodeData.normalIndex];
		}
	}

	// Fine-to-coarse down-sampling of constraints
	for(int d = maxDepth - 1; d >= (boundaryType_ == BoundaryTypeNone ? 2 : 0); --d)
		DownSample(d, sNodes_, &constraints[0]);

	// Coarse-to-fine up-sampling of coefficients
	for(int d = (boundaryType_ == BoundaryTypeNone ? 2 : 0); d < maxDepth; ++d)
		UpSample(d, sNodes_, &coefficients[0]);

	// Add the accumulated constraints from all finer depths
#pragma omp parallel for num_threads(threads_)
	for(int i = 0; i < sNodes_.nodeCount[maxDepth]; ++i)
		sNodes_.treeNodes[i]->nodeData.constraint += constraints[i];

	constraints.clear();
	shrink_to_fit(constraints);

	// Compute the contribution from all coarser depths
	for(int d = 1; d <= maxDepth; ++d) {
		DivergenceStencils stencils = SetDivergenceStencils(d, integrator, false);
		TreeNeighborKey3 neighborKey3(maxDepth);
#pragma omp parallel for num_threads(threads_) firstprivate(neighborKey3)
		for(int i = sNodes_.nodeCount[d]; i < sNodes_.nodeCount[d + 1]; ++i) {
			TreeOctNode* node = sNodes_.treeNodes[i];
			int off[3];
			node->depthAndOffset(d, off);
			Range3D range = Range3D::FullRange();
			UpdateCoarserSupportBounds(node, range);
			TreeNeighbors5 neighbors5 = neighborKey3.getNeighbors5(node->parent());

			int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 4 : 4;
			int mx = (1 << d) - mn;
			bool isInterior =
				off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
			int cx;
			int cy;
			int cz;
			Cube::FactorCornerIndex(node->parent()->childIndex(node), cx, cy, cz);
			DivergenceStencil& _stencil = stencils.at(cx, cy, cz);

			Real constraint = 0;
			for(int x = range.xStart; x != range.xEnd; ++x) {
				for(int y = range.yStart; y != range.yEnd; ++y) {
					for(int z = range.zStart; z != range.zEnd; ++z) {
						TreeOctNode* _node = neighbors5.at(x, y, z);
						if(_node && _node->nodeData.nodeIndex != -1) {
							int _d;
							int _off[3];
							_node->depthAndOffset(_d, _off);
							Point3D<Real>& normal = coefficients[_node->nodeData.nodeIndex];
							constraint += isInterior ?
								(Real)Dot(_stencil.at(x, y, z), Point3D<double>(normal)) :
								(Real)GetDivergence2(integrator, d, off, _off, true, normal);
						}
					}
				}
			}
			node->nodeData.constraint += constraint;
		}
	}

	// Set the point weights for evaluating the iso-value
#pragma omp parallel for num_threads(threads_)
	for(int i = 0; i < sNodes_.nodeCount[maxDepth + 1]; ++i) {
		TreeOctNode* temp = sNodes_.treeNodes[i];
		temp->nodeData.centerWeightContribution[OutputDensity ? 1 : 0] =
			temp->nodeData.nodeIndex < 0 || temp->nodeData.normalIndex < 0 ? 0 :
			Length(normals_[temp->nodeData.normalIndex]);
	}
	MemoryUsage();
	normals_.clear();
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::FaceEdgesFunction::operator()(TreeOctNode const* node1,
		TreeOctNode const*) const {
	if(!node1->hasChildren() && MarchingCubes::HasRoots(node1->nodeData.mcIndex)) {
		int isoTri[DIMENSION * MarchingCubes::MAX_TRIANGLES];
		int count = MarchingCubes::AddTriangleIndices(node1->nodeData.mcIndex, isoTri);

		for(int j = 0; j != count; ++j) {
			for(int k = 0; k != 3; ++k) {
				if(fIndex != Cube::FaceAdjacentToEdges(isoTri[j * 3 + k], isoTri[j * 3 + ((k + 1) % 3)]))
					continue;
				RootInfo<OutputDensity> ri1;
				RootInfo<OutputDensity> ri2;
				if(GetRootIndex(node1, isoTri[j * 3 + k], maxDepth, neighborKey3, ri1) &&
						GetRootIndex(node1, isoTri[j * 3 + ((k + 1) % 3)], maxDepth, neighborKey3, ri2)) {
					edges.push_back(std::make_pair(ri2, ri1));
					if(vertexCount.find(ri1.key) == vertexCount.end())
						vertexCount[ri1.key] = std::make_pair(ri1, 0);
					if(vertexCount.find(ri2.key) == vertexCount.end())
						vertexCount[ri2.key] = std::make_pair(ri2, 0);
					--vertexCount[ri1.key].second;
					++vertexCount[ri2.key].second;
				} else std::cerr << "Bad Edge 1: " << ri1.key << " " << ri2.key << std::endl;
			}
		}
	}
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::refineBoundary(int subdivideDepth) {
	// This implementation is somewhat tricky.
	// We would like to ensure that leaf-nodes across a subdivision boundary have the same depth.
	// We do this by calling the setNeighbors function.
	// The key is to implement this in a single pass through the leaves, ensuring that refinements
	// don't propogate. To this end, we do the minimal refinement that ensures that a cross boundary
	// neighbor, and any of its cross-boundary neighbors are all refined simultaneously.
	// For this reason, the implementation can only support nodes deeper than sDepth.
	int maxDepth = tree_.maxDepth();

	subdivideDepth = std::max(subdivideDepth, 0);
	if(boundaryType_ == BoundaryTypeNone) subdivideDepth += 2;
	subdivideDepth = std::min(subdivideDepth, maxDepth);
	int sDepth = maxDepth - subdivideDepth;
	if(boundaryType_ == BoundaryTypeNone) sDepth = std::max(2, sDepth);
	if(sDepth == 0) {
		sNodes_.set(tree_);
		return sDepth;
	}

	// Ensure that face adjacent neighbors across the subdivision boundary exist to allow for
	// a consistent definition of the iso-surface
	TreeNeighborKey3 nKey(maxDepth);
	for(TreeOctNode* leaf = tree_.nextLeaf(); leaf; leaf = tree_.nextLeaf(leaf)) {
		if(leaf->depth() <= sDepth) continue;
		int d;
		int off[3];
		leaf->depthAndOffset(d, off);
		int res = (1 << d) - 1;
		int _res = (1 << (d - sDepth)) - 1;
		int _off[3] = { off[0] & _res, off[1] & _res, off[2] & _res };
		bool boundary[3][2] = {
			{ off[0] != 0 && _off[0] == 0, off[0] != res && _off[0] == _res },
			{ off[1] != 0 && _off[1] == 0, off[1] != res && _off[1] == _res },
			{ off[2] != 0 && _off[2] == 0, off[2] != res && _off[2] == _res },
		};

		if(!boundary[0][0] && !boundary[0][1] && !boundary[1][0] && !boundary[1][1] &&
				!boundary[2][0] && !boundary[2][1]) continue;
		TreeNeighbors3& neighbors = nKey.getNeighbors3(leaf);
		int x = boundary[0][0] && !neighbors.at(0, 1, 1) ? -1 :
			boundary[0][1] && !neighbors.at(2, 1, 1) ? 1 : 0;
		int y = boundary[1][0] && !neighbors.at(1, 0, 1) ? -1 :
			boundary[1][1] && !neighbors.at(1, 2, 1) ? 1 : 0;
		int z = boundary[2][0] && !neighbors.at(1, 1, 0) ? -1 :
			boundary[2][1] && !neighbors.at(1, 1, 2) ? 1 : 0;

		if(!x && !y && !z) continue;
		bool flags[3][3][3] = {};
		// Corner case
		if(x && y && z) flags[1 + x][1 + y][1 + z] = true;
		// Edge cases
		if(x && y) flags[1 + x][1 + y][1] = true;
		if(x && z) flags[1 + x][1][1 + z] = true;
		if(y && z) flags[1][1 + y][1 + 1] = true; // TODO: maybe 1 + z?
		// Face cases
		if(x) flags[1 + x][1][1] = true;
		if(y) flags[1][1 + y][1] = true;
		if(z) flags[1][1][1 + z] = true;
		nKey.setNeighbors(leaf, flags);
	}
	sNodes_.set(tree_);
	MemoryUsage();
	return sDepth;
}

template<int Degree, bool OutputDensity>
template<class Vertex>
void Octree<Degree, OutputDensity>::GetMCIsoTriangles(Real isoValue, int subdivideDepth,
		CoredFileMeshData<Vertex>* mesh, int nonLinearFit, bool addBarycenter, bool polygonMesh) {
	CornerEvaluator2 evaluator;
	fData_.setCornerEvaluator(evaluator, 0, postDerivativeSmooth_);
	// Ensure that the subtrees are self-contained
	int sDepth = refineBoundary(subdivideDepth);

	std::vector<Vertex>* interiorVertices;
	int maxDepth = tree_.maxDepth();

	std::vector<Real> metSolution(sNodes_.nodeCount[maxDepth], 0);
#pragma omp parallel for num_threads(threads_)
	for(int i = sNodes_.nodeCount[minDepth_]; i < sNodes_.nodeCount[maxDepth]; ++i)
		metSolution[i] = sNodes_.treeNodes[i]->nodeData.solution;
	for(int d = minDepth_; d < maxDepth; ++d) UpSample(d, sNodes_, &metSolution[0]);

	// Clear the marching cube indices
#pragma omp parallel for num_threads( threads_ )
	for(int i = 0; i < sNodes_.nodeCount[maxDepth + 1]; ++i)
		sNodes_.treeNodes[i]->nodeData.mcIndex = 0;

	int offSet = 0;

	int maxCCount = sNodes_.getMaxCornerCount(sDepth, maxDepth, threads_);
	int maxECount = sNodes_.getMaxEdgeCount(&tree_, sDepth, threads_);

	RootData<OutputDensity> rootData;
	rootData.cornerValues.resize(maxCCount);
	rootData.cornerNormals.resize(maxCCount);
	rootData.interiorRoots.resize(maxECount);
	rootData.cornerValuesSet.resize(maxCCount);
	rootData.cornerNormalsSet.resize(maxCCount);
	rootData.edgesSet.resize(maxECount);
	RootData<OutputDensity> coarseRootData;
	sNodes_.setCornerTable(coarseRootData, nullptr, sDepth, threads_);
	coarseRootData.cornerValues.resize(coarseRootData.cCount());
	coarseRootData.cornerNormals.resize(coarseRootData.cCount());
	coarseRootData.cornerValuesSet.assign(coarseRootData.cCount(), 0);
	coarseRootData.cornerNormalsSet.assign(coarseRootData.cCount(), 0);
	MemoryUsage();

	TreeConstNeighborKey3 nKey(maxDepth);
	std::vector<CornerValueStencil> vStencils(maxDepth + 1);
	std::vector<CornerNormalStencil> nStencils(maxDepth + 1);
	for(int d = minDepth_; d <= maxDepth; ++d) {
		vStencils[d].stencil = SetCornerEvaluationStencil(evaluator, d);
		vStencils[d].stencils = SetCornerEvaluationStencils(evaluator, d);
		nStencils[d].stencil = SetCornerNormalEvaluationStencil(evaluator, d);
		nStencils[d].stencils = SetCornerNormalEvaluationStencils(evaluator, d);
	}

	// First process all leaf nodes at depths strictly finer than sDepth, one subtree at a time.
	for(int i = sNodes_.nodeCount[sDepth]; i != sNodes_.nodeCount[sDepth + 1]; ++i) {
		if(!sNodes_.treeNodes[i]->hasChildren()) continue;

		sNodes_.setCornerTable(rootData, sNodes_.treeNodes[i], threads_);
		sNodes_.setEdgeTable(rootData, sNodes_.treeNodes[i], threads_);
		rootData.cornerValuesSet.assign(rootData.cCount(), 0);
		rootData.cornerNormalsSet.assign(rootData.cCount(), 0);
		rootData.edgesSet.assign(rootData.eCount(), 0);
		interiorVertices = new std::vector<Vertex>();
		for(int d = maxDepth; d > sDepth; --d) {
			std::vector<TreeOctNode*> leafNodes;
			for(TreeOctNode* node = sNodes_.treeNodes[i]->nextLeaf(); node;
					node = sNodes_.treeNodes[i]->nextLeaf(node))
				if(node->depth() == d && node->nodeData.nodeIndex != -1)
					leafNodes.push_back(node);
			size_t leafNodeCount = leafNodes.size();

			// First set the corner values and associated marching-cube indices
#pragma omp parallel for num_threads(threads_) firstprivate(nKey)
			for(int j = 0; j < (int)leafNodeCount; ++j) {
				TreeOctNode* leaf = leafNodes[j];
				SetIsoCorners(isoValue, leaf, rootData, &rootData.cornerValuesSet[0],
						&rootData.cornerValues[0], nKey, metSolution, evaluator, vStencils[d].stencil,
						vStencils[d].stencils);

				// If this node shares a vertex with a coarser node, set the vertex value
				int d;
				int off[3];
				leaf->depthAndOffset(d, off);
				int res = 1 << (d - sDepth);
				off[0] %= res;
				off[1] %= res;
				off[2] %= res;
				--res;
				if(!(off[0] % res) && !(off[1] % res) && !(off[2] % res)) {
					TreeOctNode const* temp = leaf;
					while(temp->depth() != sDepth) temp = temp->parent();
					int x = off[0] == 0 ? 0 : 1;
					int y = off[1] == 0 ? 0 : 1;
					int z = off[2] == 0 ? 0 : 1;
					int c = Cube::CornerIndex(x, y, z);
					int idx = coarseRootData.cornerIndices(temp, c);
					coarseRootData.cornerValues[idx] = rootData.cornerValues[rootData.cornerIndices(leaf, c)];
					coarseRootData.cornerValuesSet[idx] = true;
				}

				// Compute the iso-vertices
				//
				if(boundaryType_ != BoundaryTypeNone || IsInset(leaf))
					SetMCRootPositions(leaf, sDepth, isoValue, nKey, rootData, interiorVertices, mesh,
							metSolution, evaluator, nStencils[d].stencil, nStencils[d].stencils,
							nonLinearFit);
			}
			// Note that this should be broken off for multi-threading as
			// the SetMCRootPositions writes to interiorPoints (with locking)
			// while GetMCIsoTriangles reads from interiorPoints (without locking)
			std::vector<Vertex> barycenters;
#pragma omp parallel for num_threads(threads_) firstprivate(nKey)
			for(int i = 0; i < (int)leafNodeCount; ++i) {
				TreeOctNode* leaf = leafNodes[i];
				if(boundaryType_ != BoundaryTypeNone || IsInset(leaf))
					GetMCIsoTriangles(leaf, nKey, mesh, rootData, interiorVertices, offSet, sDepth,
							polygonMesh, addBarycenter ? &barycenters : nullptr);
			}
			for(size_t i = 0; i != barycenters.size(); ++i) interiorVertices->push_back(barycenters[i]);
		}
		offSet = mesh->outOfCorePointCount();
		delete interiorVertices;
	}

	MemoryUsage();
	rootData.cornerNormalsSet.clear();
	shrink_to_fit(rootData.cornerNormalsSet);
	rootData.cornerValues.clear();
	shrink_to_fit(rootData.cornerValues);
	rootData.edgesSet.clear();
	shrink_to_fit(rootData.edgesSet);
	rootData.cornerValuesSet.clear();
	shrink_to_fit(rootData.cornerValuesSet);
	rootData.interiorRoots.clear();
	shrink_to_fit(rootData.interiorRoots);
	coarseRootData.interiorRoots.clear();
	coarseRootData.boundaryValues = rootData.boundaryValues;
	rootData.boundaryValues.clear();
	for(HashMap<long long, int>::iterator iter = rootData.boundaryRoots.begin();
			iter != rootData.boundaryRoots.end(); ++iter)
		coarseRootData.boundaryRoots[iter->first] = iter->second;

	for(int d = sDepth; d >= 0; --d) {
		std::vector<Vertex> barycenters;
		for(int i = sNodes_.nodeCount[d]; i != sNodes_.nodeCount[d + 1]; ++i) {
			TreeOctNode* leaf = sNodes_.treeNodes[i];
			if(leaf->hasChildren()) continue;

			// First set the corner values and associated marching-cube indices
			SetIsoCorners(isoValue, leaf, coarseRootData, &coarseRootData.cornerValuesSet[0],
					&coarseRootData.cornerValues[0], nKey, metSolution, evaluator, vStencils[d].stencil,
					vStencils[d].stencils);

			// Now compute the iso-vertices
			if(boundaryType_ != BoundaryTypeNone || IsInset(leaf)) {
				SetMCRootPositions<Vertex>(leaf, 0, isoValue, nKey, coarseRootData, nullptr, mesh,
						metSolution, evaluator, nStencils[d].stencil, nStencils[d].stencils, nonLinearFit);
				GetMCIsoTriangles<Vertex>(leaf, nKey, mesh, coarseRootData, nullptr, 0, 0, polygonMesh,
						addBarycenter ? &barycenters : nullptr);
			}
		}
	}
	MemoryUsage();
}

template<int Degree, bool OutputDensity>
Real Octree<Degree, OutputDensity>::getCenterValue(TreeConstNeighborKey3 const& neighborKey,
		TreeOctNode const* node, std::vector<Real> const& metSolution, CenterEvaluator1 const& evaluator,
		Stencil<double, 3> const& stencil, Stencil<double, 3> const& pStencil, bool isInterior) const {
	if(node->hasChildren())
		std::cerr << "[WARNING] getCenterValue assumes leaf node" << std::endl;
	Real value = 0;

	int d;
	int off[3];
	node->depthAndOffset(d, off);

	if(isInterior) {
		for(int i = 0; i != 3; ++i) {
			for(int j = 0; j != 3; ++j) {
				for(int k = 0; k != 3; ++k) {
					TreeOctNode const* n = neighborKey.neighbors(d).at(i, j, k);
					if(n) value += n->nodeData.solution * (Real)stencil.at(i, j, k);
				}
			}
		}
		if(d > minDepth_) {
			for(int i = 0; i != 3; ++i) {
				for(int j = 0; j != 3; ++j) {
					for(int k = 0; k != 3; ++k) {
						TreeOctNode const* n = neighborKey.neighbors(d - 1).at(i, j, k);
						if(n) value += metSolution[n->nodeData.nodeIndex] * (Real)pStencil.at(i, j, k);
					}
				}
			}
		}
	} else {
		for(int i = 0; i != 3; ++i) {
			for(int j = 0; j != 3; ++j) {
				for(int k = 0; k != 3; ++k) {
					TreeOctNode const* n = neighborKey.neighbors(d).at(i, j, k);
					if(n) {
						int _d;
						int _off[3];
						n->depthAndOffset(_d, _off);
						value += n->nodeData.solution * (Real)(
							evaluator.value(d, off[0], _off[0], false, false) *
							evaluator.value(d, off[1], _off[1], false, false) *
							evaluator.value(d, off[1], _off[1], false, false)); // TODO: Maybe 2?
					}
				}
			}
		}
		if(d > minDepth_) {
			for(int i = 0; i != 3; ++i) {
				for(int j = 0; j != 3; ++j) {
					for(int k = 0; k != 3; ++k) {
						const TreeOctNode* n = neighborKey.neighbors(d-1).at(i, j, k);
						if(n) {
							int _d;
							int _off[3];
							n->depthAndOffset(_d, _off);
							value += n->nodeData.solution * (Real)( // TODO: Maybe metSolution[]?
								evaluator.value(d, off[0], _off[0], false, false) *
								evaluator.value(d, off[1], _off[1], false, false) *
								evaluator.value(d, off[1], _off[1], false, false)); // TODO: Maybe 2?
						}
					}
				}
			}
		}
	}
	return value;
}

// TODO: Looks very similar to getCenterValue.
template<int Degree, bool OutputDensity>
Real Octree<Degree, OutputDensity>::getCornerValue(TreeConstNeighborKey3 const& neighborKey3,
		TreeOctNode const* node, int corner, std::vector<Real> const& metSolution,
		CornerEvaluator2 const& evaluator, Stencil<double, 3> const& stencil,
		CornerEvaluationStencil const& stencils, bool isInterior) const {
	double value = boundaryType_ == BoundaryTypeDirichlet ? -0.5 : 0;
	int d;
	int off[3];
	node->depthAndOffset(d, off);

	int cx;
	int cy;
	int cz;
	Range3D range;
	range.xStart = range.yStart = range.zStart = 0;
	range.xEnd = range.yEnd = range.zEnd = 3;
	Cube::FactorCornerIndex(corner, cx, cy, cz);
	TreeConstNeighbors3 const& neighbors = neighborKey3.neighbors(d);
	if(cx == 0) range.xEnd = 2;
	else range.xStart = 1;
	if(cy == 0) range.yEnd = 2;
	else range.yStart = 1;
	if(cz == 0) range.zEnd = 2;
	else range.zStart = 1;
	if(isInterior) {
		for(int x = range.xStart; x != range.xEnd; ++x) {
			for(int y = range.yStart; y != range.yEnd; ++y) {
				for(int z = range.zStart; z != range.zEnd; ++z) {
					TreeOctNode const* _node = neighbors.at(x, y, z);
					if(_node) value += _node->nodeData.solution * stencil.at(x, y, z);
				}
			}
		}
	} else {
		for(int x = range.xStart; x != range.xEnd; ++x) {
			for(int y = range.yStart; y != range.yEnd; ++y) {
				for(int z = range.zStart; z != range.zEnd; ++z) {
					TreeOctNode const* _node = neighbors.at(x, y, z);
					if(_node) {
						int _d;
						int _off[3];
						_node->depthAndOffset(_d, _off);
						value += _node->nodeData.solution *
							evaluator.value(d, off[0], cx, _off[0], false, false) *
							evaluator.value(d, off[1], cy, _off[1], false, false) *
							evaluator.value(d, off[2], cz, _off[2], false, false);
					}
				}
			}
		}
	}
	if(d > minDepth_) {
		int _corner = node->parent()->childIndex(node);
		int _cx;
		int _cy;
		int _cz;
		Cube::FactorCornerIndex(_corner, _cx, _cy, _cz);
		if(cx !=_cx) {
			range.xStart = 0;
			range.xEnd = 3;
		}
		if(cy !=_cy) {
			range.yStart = 0;
			range.yEnd = 3;
		}
		if(cz !=_cz) {
			range.zStart = 0;
			range.zEnd = 3;
		}
		TreeConstNeighbors3 const& neighbors = neighborKey3.neighbors(d - 1);
		if(isInterior) {
			for(int x = range.xStart; x != range.xEnd; ++x) {
				for(int y = range.yStart; y != range.yEnd; ++y) {
					for(int z = range.zStart; z != range.zEnd; ++z) {
						TreeOctNode const* _node = neighbors.at(x, y, z);
						if(_node) value +=
							metSolution[_node->nodeData.nodeIndex] * stencils.at(_cx, _cy, _cz).at(x, y, z);
					}
				}
			}
		} else {
			for(int x = range.xStart; x != range.xEnd; ++x) {
				for(int y = range.yStart; y != range.yEnd; ++y) {
					for(int z = range.zStart; z != range.zEnd; ++z) {
						const TreeOctNode* _node = neighbors.at(x, y, z);
						if(_node) {
							int _d;
							int _off[3];
							_node->depthAndOffset(_d, _off);
							value += metSolution[_node->nodeData.nodeIndex] *
								evaluator.value(d, off[0], cx, _off[0], false, true) *
								evaluator.value(d, off[1], cy, _off[1], false, true) *
								evaluator.value(d, off[2], cz, _off[2], false, true);
						}
					}
				}
			}
		}
	}
	return value;
}

// TODO: Looks similar to getCornerValue
template<int Degree, bool OutputDensity>
Point3D<Real> Octree<Degree, OutputDensity>::getCornerNormal(TreeConstNeighbors5 const& neighbors5,
		TreeConstNeighbors5 const& pNeighbors5, TreeOctNode const* node, int corner,
		std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
		Stencil<Point3D<double>, 5> const& nStencil, CornerNormalEvaluationStencil const& nStencils,
		bool isInterior) const {
	Point3D<double> normal;

	int d;
	int off[3];
	node->depthAndOffset(d, off);

	int cx;
	int cy;
	int cz;
	Range3D range = Range3D::FullRange();
	Cube::FactorCornerIndex(corner, cx, cy, cz);
	if(cx == 0) range.xEnd = 4;
	else range.xStart = 1;
	if(cy == 0) range.yEnd = 4;
	else range.yStart = 1;
	if(cz == 0) range.zEnd = 4;
	else range.zStart = 1;
	if(isInterior) {
		for(int x = range.xStart; x != range.xEnd; ++x) {
			for(int y = range.yStart; y != range.yEnd; ++y) {
				for(int z = range.zStart; z != range.zEnd; ++z) {
					TreeOctNode const* _node = neighbors5.at(x, y, z);
					if(_node) normal += nStencil.at(x, y, z) * _node->nodeData.solution;
				}
			}
		}
	} else {
		for(int x = range.xStart; x != range.xEnd; ++x) {
			for(int y = range.yStart; y != range.yEnd; ++y) {
				for(int z = range.zStart; z != range.zEnd; ++z) {
					TreeOctNode const* _node = neighbors5.at(x, y, z);
					if(_node) {
						int _d;
						int _off[3];
						_node->depthAndOffset(_d, _off);
						double v[] = { evaluator.value(d, off[0], cx, _off[0], false, false),
							evaluator.value(d, off[1], cy, _off[1], false, false),
							evaluator.value(d, off[2], cz, _off[2], false, false) };
						double dv[] = { evaluator.value(d, off[0], cx, _off[0], true, false),
							evaluator.value(d, off[1], cy, _off[1], true, false),
							evaluator.value(d, off[2], cz, _off[2], true, false) };
						normal +=
							Point3D<double>(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]) *
							_node->nodeData.solution;
					}
				}
			}
		}
	}
	if(d > minDepth_) {
		int _cx;
		int _cy;
		int _cz;
		int _corner = node->parent()->childIndex(node);
		Cube::FactorCornerIndex(_corner, _cx, _cy, _cz);
		if(cx !=_cx) {
			range.xStart = 0;
			range.xEnd = 5;
		}
		if(cy !=_cy) {
			range.yStart = 0;
			range.yEnd = 5;
		}
		if(cz !=_cz) {
			range.zStart = 0;
			range.zEnd = 5;
		}
		if(isInterior) {
			for(int x = range.xStart; x != range.xEnd; ++x) {
				for(int y = range.yStart; y != range.yEnd; ++y) {
					for(int z = range.zStart; z != range.zEnd; ++z) {
						TreeOctNode const* _node = pNeighbors5.at(x, y, z);
						if(_node) normal += nStencils.at(_cx, _cy, _cz).at(x, y, z) *
							metSolution[_node->nodeData.nodeIndex];
					}
				}
			}
		} else {
			for(int x = range.xStart; x != range.xEnd; ++x) {
				for(int y = range.yStart; y != range.yEnd; ++y) {
					for(int z = range.zStart; z != range.zEnd; ++z) {
						TreeOctNode const* _node = pNeighbors5.at(x, y, z);
						if(_node) {
							int _d;
							int _off[3];
							_node->depthAndOffset(_d, _off);
							double v[] = { evaluator.value(d, off[0], cx, _off[0], false, true),
								evaluator.value(d, off[1], cy, _off[1], false, true),
								evaluator.value(d, off[2], cz, _off[2], false, true) };
							double dv[] = { evaluator.value(d, off[0], cx, _off[0], true, true),
								evaluator.value(d, off[1], cy, _off[1], true, true),
								evaluator.value(d, off[2], cz, _off[2], true, true) };
							normal +=
								Point3D<double>(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]) *
								metSolution[_node->nodeData.nodeIndex];
						}
					}
				}
			}
		}
	}
	return Point3D<Real>(normal);
}

template<int Degree, bool OutputDensity>
Real Octree<Degree, OutputDensity>::GetIsoValue() const {
	Real isoValue = 0;
	Real weightSum = 0;
	int maxDepth = tree_.maxDepth();

	CenterEvaluator1 evaluator;
	fData_.setCenterEvaluator(evaluator, 0, 0);
	std::vector<CenterValueStencil> vStencils(maxDepth + 1);
	for(int d = minDepth_; d <= maxDepth; ++d) {
		vStencils[d].stencil = SetCenterEvaluationStencil(evaluator, d);
		vStencils[d].stencils = SetCenterEvaluationStencils(evaluator, d);
	}
	std::vector<Real> metSolution(sNodes_.nodeCount[maxDepth], 0);
	std::vector<Real> centerValues(sNodes_.nodeCount[maxDepth + 1]);
#pragma omp parallel for num_threads(threads_)
	for(int i = sNodes_.nodeCount[minDepth_]; i < sNodes_.nodeCount[maxDepth]; ++i)
		metSolution[i] = sNodes_.treeNodes[i]->nodeData.solution;
	for(int d = minDepth_; d < maxDepth; ++d)
		UpSample(d, sNodes_, &metSolution[0] + sNodes_.nodeCount[d - 1],
				&metSolution[0] + sNodes_.nodeCount[d]);
	for(int d = maxDepth; d >= minDepth_; --d) {
		TreeConstNeighborKey3 nKey(d);
#pragma omp parallel for num_threads(threads_) reduction(+ : isoValue, weightSum) firstprivate(nKey)
		for(int i = sNodes_.nodeCount[d]; i < sNodes_.nodeCount[d + 1]; ++i) {
			TreeOctNode* node = sNodes_.treeNodes[i];
			Real value = 0;
			if(node->hasChildren()) {
				for(unsigned c = 0; c < Cube::CORNERS; ++c)
					value += centerValues[node->child(c)->nodeData.nodeIndex];
				value /= Cube::CORNERS;
			} else {
				nKey.getNeighbors3(node);
				int c = 0;
				int x;
				int y;
				int z;
				if(node->parent()) c = node->parent()->childIndex(node);
				Cube::FactorCornerIndex(c, x, y, z);

				int d;
				int off[3];
				node->depthAndOffset(d, off);
				int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 2 : 2;
				int mx = (1 << d) - mn;
				bool isInterior =
					off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;

				value = getCenterValue(nKey, node, metSolution, evaluator, vStencils[d].stencil,
						vStencils[d].stencils.at(x, y, z), isInterior);
			}
			centerValues[i] = value;
			Real w = node->nodeData.centerWeightContribution[OutputDensity ? 1 : 0];
			if(w != 0) {
				isoValue += value * w;
				weightSum += w;
			}
		}
	}
	Real r = boundaryType_ == BoundaryTypeDirichlet ? 0.5 : 0;
	return isoValue / weightSum - r;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::SetIsoCorners(Real isoValue, TreeOctNode* leaf,
		CornerTableData& cData, char* valuesSet, Real* values, TreeConstNeighborKey3& nKey,
		std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
		CornerEvaluationStencil const& stencil, CornerEvaluationStencils const& stencils) {
	Real cornerValues[Cube::CORNERS];
	typename SortedTreeNodes< OutputDensity >::CornerIndices const& cIndices = cData[leaf];

	int d;
	int off[3];
	leaf->depthAndOffset(d, off);
	int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 2 : 2;
	int mx = (1 << d) - mn;
	bool isInterior =
		off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
	nKey.getNeighbors3(leaf);
	for(unsigned c = 0; c != Cube::CORNERS; ++c) {
		int vIndex = cIndices[c];
		if(valuesSet[vIndex]) cornerValues[c] = values[vIndex];
		else {
			int x;
			int y;
			int z;
			Cube::FactorCornerIndex(c, x, y, z);
			cornerValues[c] = getCornerValue(nKey, leaf, c, metSolution, evaluator,
					stencil.at(x, y, z), stencils.at(x, y, z), isInterior);
			values[vIndex] = cornerValues[c];
			valuesSet[vIndex] = 1;
		}
	}
	leaf->nodeData.mcIndex = MarchingCubes::GetIndex(cornerValues, isoValue);

	// Set the marching cube indices for all interior nodes.
	if(leaf->parent()) {
		TreeOctNode* parent = leaf->parent();
		int c = leaf->parent()->childIndex(leaf);
		int mcid = leaf->nodeData.mcIndex & (1 << MarchingCubes::cornerMap[c]);

		if(mcid) {
#ifdef WIN32
			InterlockedOr( (volatile unsigned long long*)&(parent->nodeData.mcIndex) , mcid );
#else
#pragma omp atomic
			parent->nodeData.mcIndex |= mcid;
#endif
			while(1) {
				if(parent->parent() && parent->parent()->depth() >= minDepth_ &&
						parent->parent()->childIndex(parent) == c) {
#ifdef WIN32
					InterlockedOr( (volatile unsigned long long*)&(parent->parent()->nodeData.mcIndex) , mcid );
#else
#pragma omp atomic
					parent->parent()->nodeData.mcIndex |= mcid;
#endif
					parent = parent->parent();
				}
				else break;
			}
		}
	}
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::IsBoundaryFace(TreeOctNode const* node, int faceIndex, int sDepth) {
	if(sDepth < 0) return 0;
	if(node->depth() <= sDepth) return 1;

	int dir;
	int offset;
	Cube::FactorFaceIndex(faceIndex, dir, offset);
	int d;
	int o[3];
	node->depthAndOffset(d, o);

	int idx = (o[dir] << 1) + (offset << 1);
	return !(idx % (2 << (d - sDepth)));
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::IsBoundaryEdge(TreeOctNode const* node, int edgeIndex, int sDepth) {
	int dir;
	int x;
	int y;
	Cube::FactorEdgeIndex(edgeIndex, dir, x, y);
	return IsBoundaryEdge(node, dir, x, y, sDepth);
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::IsBoundaryEdge(TreeOctNode const* node, int dir, int x, int y,
		int sDepth) {
	if(sDepth < 0 ) return 0;
	if(node->depth() <= sDepth) return 1;

	int d;
	int o[3];
	node->depthAndOffset( d , o );

	int idx1 = 0;
	int idx2 = 0;
	switch(dir) {
		case 0:
			idx1 = o[1] + x;
			idx2 = o[2] + y;
			break;
		case 1:
			idx1 = o[0] + x;
			idx2 = o[2] + y;
			break;
		case 2:
			idx1 = o[0] + x;
			idx2 = o[1] + y;
			break;
	}
	int mask = 1 << (d - sDepth);
	return !(idx1 % mask) || !(idx2 % mask);
}

//////////////////////////////////////////////////////////////////////////////////////
// The assumption made when calling this code is that the edge has at most one root //
//////////////////////////////////////////////////////////////////////////////////////

// TODO: Required by GetRoot
void SetVertexValue(PlyVertex<Real>&, Real) { }
void SetVertexValue(PlyValueVertex<Real>& vertex, Real value) { vertex.value = value; }

template<int Degree, bool OutputDensity>
template<class Vertex>
int Octree<Degree, OutputDensity>::GetRoot(RootInfo<OutputDensity> const& ri, Real isoValue,
		TreeConstNeighborKey3& neighborKey3, Vertex& vertex, RootData<OutputDensity>& rootData,
		int sDepth, std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
		CornerNormalEvaluationStencil const& nStencil, CornerNormalEvaluationStencils const& nStencils,
		bool nonLinearFit) {
	if(!MarchingCubes::HasRoots(ri.node->nodeData.mcIndex)) return 0;
	if(!MarchingCubes::HasEdgeRoots(ri.node->nodeData.mcIndex, ri.edgeIndex)) return 0;

	int c1;
	int c2;
	Cube::EdgeCorners(ri.edgeIndex, c1, c2);

	int o;
	int i1;
	int i2;
	Cube::FactorEdgeIndex(ri.edgeIndex, o, i1, i2);

	long long key1 = VertexData<OutputDensity>::CornerIndex(ri.node, c1, fData_.depth());
	long long key2 = VertexData<OutputDensity>::CornerIndex(ri.node, c2, fData_.depth());

	bool isBoundary = IsBoundaryEdge(ri.node, ri.edgeIndex, sDepth) != 0;
	bool haveKey1;
	bool haveKey2;
	std::pair<Real, Point3D<Real> > keyValue1;
	std::pair<Real, Point3D<Real> > keyValue2;
	int iter1 = rootData.cornerIndices(ri.node, c1);
	int iter2 = rootData.cornerIndices(ri.node, c2);
	keyValue1.first = rootData.cornerValues[iter1];
	keyValue2.first = rootData.cornerValues[iter2];
	if(isBoundary) {
#pragma omp critical (normal_hash_access)
		{
			haveKey1 = rootData.boundaryValues.find(key1) != rootData.boundaryValues.end();
			haveKey2 = rootData.boundaryValues.find(key2) != rootData.boundaryValues.end();
			if(haveKey1) keyValue1 = rootData.boundaryValues[key1];
			if(haveKey2) keyValue2 = rootData.boundaryValues[key2];
		}
	} else {
		haveKey1 = rootData.cornerNormalsSet[iter1] != 0;
		haveKey2 = rootData.cornerNormalsSet[iter2] != 0;
		if(haveKey1) keyValue1.second = rootData.cornerNormals[iter1];
		if(haveKey2) keyValue2.second = rootData.cornerNormals[iter2];
	}
	TreeConstNeighbors5 neighbors5;
	TreeConstNeighbors5 pNeighbors5;
	bool isInterior = false;
	if(!haveKey1 || !haveKey2) {
		neighbors5 = neighborKey3.getNeighbors5(ri.node);
		if(ri.node->parent()) pNeighbors5 = neighborKey3.getNeighbors5(ri.node->parent());
		int d;
		int off[3];
		ri.node->depthAndOffset(d, off);
		int mn = boundaryType_ == BoundaryTypeNone ? (1 << (d - 2)) + 2 : 2;
		int mx = (1 << d) - mn;
		isInterior =
			off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx;
	}
	int c1x;
	int c1y;
	int c1z;
	Cube::FactorCornerIndex(c1, c1x, c1y, c1z);
	int c2x;
	int c2y;
	int c2z;
	Cube::FactorCornerIndex(c2, c2x, c2y, c2z);
	if(!haveKey1)
		keyValue1.second = getCornerNormal(neighbors5, pNeighbors5, ri.node, c1, metSolution, evaluator,
				nStencil.at(c1x, c1y, c1z), nStencils.at(c1x, c1y, c1z), isInterior);
	if(!haveKey2)
		keyValue2.second = getCornerNormal(neighbors5, pNeighbors5, ri.node, c2, metSolution, evaluator,
				nStencil.at(c2x, c2y, c2z), nStencils.at(c2x, c2y, c2z), isInterior);
	Point3D<Real> n[2] = { keyValue1.second, keyValue2.second };
	double x0 = keyValue1.first;
	double x1 = keyValue2.first;

	if(!haveKey1 || !haveKey2) {
		if(isBoundary) {
#pragma omp critical (normal_hash_access)
			{
				if(!haveKey1) rootData.boundaryValues[key1] = keyValue1;
				if(!haveKey2) rootData.boundaryValues[key2] = keyValue2;
			}
		} else {
			if(!haveKey1) {
				rootData.cornerNormals[iter1] = keyValue1.second;
				rootData.cornerNormalsSet[iter1] = 1;
			}
			if(!haveKey2) {
				rootData.cornerNormals[iter2] = keyValue2.second;
				rootData.cornerNormalsSet[iter2] = 1;
			}
		}
	}

	Point3D<Real> c;
	Real width;
	ri.node->centerAndWidth(c, width);
	Real center = c[o];
	for(int i = 0; i != DIMENSION; ++i) {
		n[0][i] *= width;
		n[1][i] *= width;
	}

	Point3D<Real> position;
	switch(o) {
		case 0:
			position[1] = c[1] - width / 2 + width * i1;
			position[2] = c[2] - width / 2 + width * i2;
			break;
		case 1:
			position[0] = c[0] - width / 2 + width * i1;
			position[2] = c[2] - width / 2 + width * i2;
			break;
		case 2:
			position[0] = c[0] - width / 2 + width * i1;
			position[1] = c[1] - width / 2 + width * i2;
			break;
	}
	double dx0 = n[0][o];
	double dx1 = n[1][o];

	// The scaling will turn the Hermite Spline into a quadratic
	double scl = (x1 - x0) / ((dx1 + dx0) / 2);
	dx0 *= scl;
	dx1 *= scl;

	// Hermite Spline
	double coefficients[] = { x0, dx0, 3 * (x1 - x0) - dx1 - 2 * dx0 };
	Polynomial<2> P(coefficients);

	std::vector<double> roots = P.getSolutions(isoValue, EPSILON);
	int rCount = 0;
	Real averageRoot = 0;
	for(int i = 0; i != (int)roots.size(); ++i) {
		if(roots[i] >= 0 && roots[i] <= 1) {
			averageRoot += (Real)roots[i];
			++rCount;
		}
	}
	if(rCount && nonLinearFit) averageRoot /= rCount;
	else averageRoot = (Real)((x0 - isoValue) / (x0 - x1));
	if(averageRoot < 0 || averageRoot > 1) {
		std::cerr << "[WARNING] Bad average root: " << averageRoot << std::endl << 
			"\t(" << x0 << " " << x1 << ") , (" << dx0 << " " << dx1 << ") (" <<
			isoValue << ")" << std::endl;
		averageRoot = clamp(averageRoot, 0, 1);
	}
	position[o] = center - width / 2 + width * averageRoot;
	vertex.point = position;
	if(OutputDensity) {
		Real depth;
		Real weight;
		TreeOctNode const* temp = ri.node;
		while(temp->depth() > splatDepth_) temp = temp->parent();
		GetSampleDepthAndWeight(temp, position, GetRootGetNeighborsFunction(neighborKey3),
				samplesPerNode_, depth, weight);
		SetVertexValue(vertex, depth);
	}
	return 1;
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::GetRootIndex(TreeOctNode const* node, int edgeIndex, int maxDepth,
		TreeConstNeighborKey3& neighborKey3, RootInfo<OutputDensity>& ri) {
	if(node->nodeData.nodeIndex == -1)
		std::cerr << "[WARNING] Called GetRootIndex with bad node" << std::endl;
	// The assumption is that the super-edge has a root along it. 
	if(!(MarchingCubes::edgeMask[node->nodeData.mcIndex] & (1 << edgeIndex))) return 0;

	int f1;
	int f2;
	Cube::FacesAdjacentToEdge(edgeIndex, f1, f2);

	TreeOctNode const* finest = node;
	int finestIndex = edgeIndex;
	if(node->depth() < maxDepth && !node->hasChildren()) {
		TreeConstNeighbors3& neighbors = neighborKey3.getNeighbors3(node);
		int x;
		int y;
		int z;
		Cube::FactorFaceIndexXYZ(f1, x, y, z);
		TreeOctNode const* temp = neighbors.at(x + 1, y + 1, z + 1);
		if(temp && temp->nodeData.nodeIndex != -1 && temp->hasChildren()) {
			finest = temp;
			finestIndex = Cube::FaceReflectEdgeIndex(edgeIndex, f1);
		} else {
			int x;
			int y;
			int z;
			Cube::FactorFaceIndexXYZ(f2, x, y, z);
			temp = neighbors.at(x + 1, y + 1, z + 1);
			if(temp && temp->nodeData.nodeIndex != -1 && temp->hasChildren()) {
				finest = temp;
				finestIndex = Cube::FaceReflectEdgeIndex(edgeIndex, f2);
			} else {
				int orientation;
				int d1;
				int d2;
				Cube::FactorEdgeIndex(edgeIndex, orientation, d1, d2);
				switch(orientation) {
					case 0: temp = neighbors.at(1, d1 << 1, d2 << 1); break;
					case 1: temp = neighbors.at(d1 << 1, 1, d2 << 1); break;
					case 2: temp = neighbors.at(d1 << 1, d2 << 1, 1); break;
				}
				if(temp && temp->nodeData.nodeIndex != -1 && temp->hasChildren()) {
					finest = temp;
					finestIndex = Cube::EdgeReflectEdgeIndex(edgeIndex);
				}
			}
		}
	}

	int c1;
	int c2;
	Cube::EdgeCorners(finestIndex, c1, c2);
	if(finest->hasChildren()) {
		if(GetRootIndex(finest->child(c1), finestIndex, maxDepth, neighborKey3, ri)) return 1;
		else if(GetRootIndex(finest->child(c2), finestIndex, maxDepth, neighborKey3, ri)) return 1;
		else {
			int d1;
			int off1[3];
			node->depthAndOffset(d1, off1);
			int d2;
			int off2[3];
			finest->depthAndOffset(d2, off2);
			std::cerr << "[WARNING] Couldn't find root index with either child [" << d1 <<
				"] (" << off1[0] << " " << off1[1] << " " << off1[2] << ") -> [" << d2 <<
				"] (" << off2[0] << " " << off2[1] << " " << off2[2] << ") (" <<
				node->hasChildren() << " " << finest->hasChildren() << ")" << std::endl;
			std::cout << "\t";
			for(int i = 0; i != 8; ++i) {
				if(node->nodeData.mcIndex & (1 << i)) std::cout << "1";
				else std::cout << "0";
			}
			std::cout << "\t";
			for(int i = 0; i != 8; ++i) {
				if(finest->nodeData.mcIndex & (1 << i)) std::cout << "1";
				else std::cout << "0";
			}
			std::cout << std::endl;
			return 0;
		}
	} else {
		int o;
		int i1;
		int i2;
		Cube::FactorEdgeIndex(finestIndex, o, i1, i2);
		int d;
		int off[3];
		finest->depthAndOffset(d, off);
		ri.node = finest;
		ri.edgeIndex = finestIndex;
		int eIndex[2] = {};
		int offset = BinaryNode<Real>::CenterIndex(d, off[o]);
		switch(o) {
			case 0:
				eIndex[0] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, off[1], i1);
				eIndex[1] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, off[2], i2);
				break;
			case 1:
				eIndex[0] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, off[0], i1);
				eIndex[1] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, off[2], i2);
				break;
			case 2:
				eIndex[0] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, off[0], i1);
				eIndex[1] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, off[1], i2);
				break;
		}
		ri.key = (long long)(o) | (long long)(eIndex[0]) << 5 |
			(long long)(eIndex[1]) << 25 | (long long)(offset) << 45;
		return 1;
	}
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::GetRootPair(RootInfo<OutputDensity> const& ri, int maxDepth,
		TreeConstNeighborKey3& nKey3, RootInfo<OutputDensity>& pair) {
	int c1;
	int c2;
	Cube::EdgeCorners(ri.edgeIndex, c1, c2);
	for(TreeOctNode const* node = ri.node; node->parent(); node = node->parent()) {
		int c = node->parent()->childIndex(node);
		if(c != c1 && c != c2) return 0;
		if(!MarchingCubes::HasEdgeRoots(node->parent()->nodeData.mcIndex, ri.edgeIndex))
			return GetRootIndex(node->parent()->child(c == c1 ? c2 : c1), ri.edgeIndex, maxDepth, nKey3, pair);
	}
	return 0;
}

template<int Degree, bool OutputDensity>
int Octree<Degree, OutputDensity>::GetRootIndex(RootInfo<OutputDensity> const& ri,
		RootData<OutputDensity>& rootData, CoredPointIndex& index) {
	HashMap<long long, int>::iterator rootIter = rootData.boundaryRoots.find(ri.key);
	if(rootIter != rootData.boundaryRoots.end()) {
		index.inCore = 1;
		index.index = rootIter->second;
		return 1;
	} else if(!rootData.interiorRoots.empty()) {
		int eIndex = rootData.edgeIndices(ri.node, ri.edgeIndex);
		if(rootData.edgesSet[eIndex]) {
			index.inCore = 0;
			index.index = rootData.interiorRoots[eIndex];
			return 1;
		}
	}
	return 0;
}

template<int Degree, bool OutputDensity>
template<class Vertex>
int Octree<Degree, OutputDensity>::SetMCRootPositions(TreeOctNode* node, int sDepth, Real isoValue,
		TreeConstNeighborKey3& neighborKey3, RootData<OutputDensity>& rootData,
		std::vector<Vertex>* interiorVertices, CoredFileMeshData<Vertex>* mesh,
		std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
		CornerNormalEvaluationStencil const& nStencil, CornerNormalEvaluationStencils const& nStencils,
		bool nonLinearFit) {
	int count = 0;
	if(!MarchingCubes::HasRoots(node->nodeData.mcIndex)) return 0;
	for(int i = 0; i != DIMENSION; ++i) {
		for(int j = 0; j != 2; ++j) {
			for(int k = 0; k != 2; ++k) {
				int eIndex = Cube::EdgeIndex(i, j, k);
				RootInfo<OutputDensity> ri;
				Vertex vertex;
				if(!GetRootIndex(node, eIndex, fData_.depth(), neighborKey3, ri)) continue;
				if(rootData.interiorRoots.empty() || IsBoundaryEdge(node, i, j, k, sDepth)) {
					HashMap<long long, int>::iterator iter;
					HashMap<long long, int>::iterator end;
					// Check if the root has already been set
#pragma omp critical (boundary_roots_hash_access)
					{
						iter = rootData.boundaryRoots.find(ri.key);
						end  = rootData.boundaryRoots.end();
					}
					if(iter != end) continue;
					// Get the root information
					GetRoot(ri, isoValue, neighborKey3, vertex, rootData, sDepth, metSolution, evaluator,
							nStencil, nStencils, nonLinearFit);
					vertex.point = vertex.point * scale_ + center_;
					// Add the root if it hasn't been added already
#pragma omp critical (boundary_roots_hash_access)
					{
						iter = rootData.boundaryRoots.find(ri.key);
						end = rootData.boundaryRoots.end();
						if(iter == end) {
							mesh->addInCorePoint(vertex);
							rootData.boundaryRoots[ri.key] = (int)mesh->inCorePointCount() - 1;
						}
					}
					if(iter == end) ++count;
				} else {
					int nodeEdgeIndex = rootData.edgeIndices(ri.node, ri.edgeIndex);
					if(rootData.edgesSet[nodeEdgeIndex]) continue;
					// Get the root information
					GetRoot(ri, isoValue, neighborKey3, vertex, rootData, sDepth, metSolution, evaluator,
							nStencil, nStencils, nonLinearFit);
					vertex.point = vertex.point * scale_ + center_;
					// Add the root if it hasn't been added already
#pragma omp critical (add_point_access)
					{
						if(!rootData.edgesSet[nodeEdgeIndex]) {
							rootData.interiorRoots[nodeEdgeIndex] = mesh->addOutOfCorePoint(vertex);
							interiorVertices->push_back(vertex);
							rootData.edgesSet[nodeEdgeIndex] = 1;
							++count;
						}
					}
				}
			}
		}
	}
	return count;
}

template<int Degree, bool OutputDensity>
void Octree<Degree, OutputDensity>::GetMCIsoEdges(TreeOctNode* node, TreeConstNeighborKey3& neighborKey3,
		int sDepth, edges_t& edges) {
	vertex_count_t vertexCount;

	FaceEdgesFunction fef(fData_.depth(), edges, vertexCount, neighborKey3);
	TreeOctNode const* _temp[Cube::NEIGHBORS];
	TreeConstNeighbors3& neighbors = neighborKey3.getNeighbors3(node);
	for(unsigned f = 0; f != Cube::NEIGHBORS; ++f) {
		int x;
		int y;
		int z;
		Cube::FactorFaceIndexXYZ(f, x, y, z);
		_temp[f] = neighbors.at(x + 1, y + 1, z + 1);
	}
	int isoTri[DIMENSION * MarchingCubes::MAX_TRIANGLES];
	int count = MarchingCubes::AddTriangleIndices(node->nodeData.mcIndex, isoTri);
	for(int fIndex = 0; fIndex < (int)Cube::NEIGHBORS; ++fIndex) {
		int ref = Cube::FaceReflectFaceIndex(fIndex, fIndex);
		fef.setFIndex(ref);
		TreeOctNode const* temp = _temp[fIndex];
		// If the face neighbor exists and has higher resolution than the current node,
		// get the iso-curve from the neighbor
		if(temp && temp->nodeData.nodeIndex != -1 && temp->hasChildren() &&
				!IsBoundaryFace(node, fIndex, sDepth))
			temp->processNodeFaces(temp, fef, ref);
		// Otherwise, get it from the node
		else {
			RootInfo<OutputDensity> ri1;
			RootInfo<OutputDensity> ri2;
			for(int j = 0; j != count; ++j) {
				for(int k = 0; k != 3; ++k) {
					int tri1 = isoTri[j * 3 + k];
					int tri2 = isoTri[j * 3 + (k + 1) % 3];
					if(fIndex != Cube::FaceAdjacentToEdges(tri1, tri2)) continue;
					if(GetRootIndex(node, tri1, fData_.depth(), neighborKey3, ri1) &&
							GetRootIndex(node, tri2, fData_.depth(), neighborKey3, ri2)) {
						edges.push_back(std::make_pair(ri1, ri2));
						if(vertexCount.find(ri1.key) == vertexCount.end()) {
							vertexCount[ri1.key].first = ri1;
							vertexCount[ri1.key].second = 0;
						}
						if(vertexCount.find(ri2.key) == vertexCount.end()) {
							vertexCount[ri2.key].first = ri2;
							vertexCount[ri2.key].second = 0;
						}
						++vertexCount[ri1.key].second;
						--vertexCount[ri2.key].second;
					} else std::cerr << "Bad Edge 2: " << ri1.key << " " << ri2.key << "\t" <<
						MarchingCubes::HasEdgeRoots(node->nodeData.mcIndex, tri1) << " " <<
						MarchingCubes::HasEdgeRoots(node->nodeData.mcIndex, tri2) << std::endl;
				}
			}
		}
	}
	for(int i = 0; i != (int)edges.size(); ++i) {
		if(vertexCount.find(edges[i].first.key) == vertexCount.end())
			std::cout << "Could not find vertex: " << edges[i].first.key << std::endl;
		else if(vertexCount[edges[i].first.key].second) {
			RootInfo<OutputDensity> ri;
			GetRootPair(vertexCount[edges[i].first.key].first, fData_.depth(), neighborKey3, ri);
			if(vertexCount.find(ri.key) == vertexCount.end()) {
				int d;
				int off[3];
				node->depthAndOffset(d, off);
				std::cout << "Vertex pair not in list 1 (" << ri.key << ") " << 
					IsBoundaryEdge(ri.node, ri.edgeIndex, sDepth) << "\t[" << d << "] (" <<
					off[0] << " " << off[1] << " " << off[2] << ")" << std::endl;
			} else {
				edges.push_back(std::make_pair(ri, edges[i].first));
				++vertexCount[ri.key].second;
				--vertexCount[edges[i].first.key].second;
			}
		}

		if(vertexCount.find(edges[i].second.key) == vertexCount.end())
			std::cerr << "Could not find vertex: " << edges[i].second.key << std::endl;
		else if(vertexCount[edges[i].second.key].second) {
			RootInfo<OutputDensity> ri;
			GetRootPair(vertexCount[edges[i].second.key].first, fData_.depth(), neighborKey3, ri);
			if(vertexCount.find(ri.key) == vertexCount.end()) {
				int d;
				int off[3];
				node->depthAndOffset(d, off);
				std::cout << "Vertex pair not in list 2\t[" << d << "] (" <<
					off[0] << " " << off[1] << " " << off[2] << ")" << std::endl;
			} else {
				edges.push_back(std::make_pair(edges[i].second, ri));
				--vertexCount[ri.key].second;
				++vertexCount[edges[i].second.key].second;
			}
		}
	}
}

template<int Degree, bool OutputDensity>
template<class Vertex>
int Octree<Degree, OutputDensity>::GetMCIsoTriangles(TreeOctNode* node,
		TreeConstNeighborKey3& neighborKey3, CoredFileMeshData<Vertex>* mesh,
		RootData<OutputDensity>& rootData, std::vector<Vertex>* interiorVertices, int offSet,
		int sDepth, bool polygonMesh, std::vector<Vertex>* barycenters) {
	edges_t edges;
	GetMCIsoEdges(node, neighborKey3, sDepth, edges);

	std::vector<edges_t> edgeLoops = GetEdgeLoops(edges);

	int tris = 0;
	for(int i = 0; i != (int)edgeLoops.size(); ++i) {
		std::vector<CoredPointIndex> edgeIndices;
		for(int j = edgeLoops[i].size() - 1; j >= 0; --j) {
			CoredPointIndex p;
			if(!GetRootIndex(edgeLoops[i][j].first, rootData, p))
				std::cout << "Bad Point Index" << std::endl;
			else edgeIndices.push_back(p);
		}
		tris += AddTriangles(mesh, edgeIndices, interiorVertices, offSet, polygonMesh, barycenters);
	}
	return tris;
}

template<int Degree, bool OutputDensity>
std::vector<typename Octree<Degree, OutputDensity>::edges_t> Octree<Degree, OutputDensity>::GetEdgeLoops(
		edges_t& edges) {
	std::vector<edges_t> loops;
	int loopSize = 0;
	while(edges.size()) {
		edges_t front;
		edges_t back;
		edge_t e = edges[0];
		loops.resize(loopSize + 1);
		edges[0] = edges.back();
		edges.pop_back();
		long long frontIdx = e.second.key;
		long long backIdx = e.first.key;
		for(int j = edges.size() - 1; j >= 0; --j) {
			if(edges[j].first.key == frontIdx || edges[j].second.key == frontIdx) {
				edge_t temp = edges[j].first.key == frontIdx ? edges[j] :
					std::make_pair(edges[j].second, edges[j].first);
				frontIdx = temp.second.key;
				front.push_back(temp);
				edges[j] = edges.back();
				edges.pop_back();
				j = edges.size();
			} else if(edges[j].first.key == backIdx || edges[j].second.key == backIdx) {
				edge_t temp = edges[j].second.key == backIdx ? edges[j] :
					std::make_pair(edges[j].second, edges[j].first);
				backIdx = temp.first.key;
				back.push_back(temp);
				edges[j] = edges.back();
				edges.pop_back();
				j = edges.size();
			}
		}
		for(int j = back.size() - 1; j >= 0; --j)
			loops[loopSize].push_back(back[j]);
		loops[loopSize].push_back(e);
		for(int j = 0; j != (int)front.size(); ++j)
			loops[loopSize].push_back(front[j]);
		loopSize++;
	}
	return loops;
}

template<int Degree, bool OutputDensity>
template<class Vertex>
int Octree<Degree, OutputDensity>::AddTriangles(CoredFileMeshData<Vertex>* mesh,
		std::vector<CoredPointIndex>& edges, std::vector<Vertex>* interiorVertices, int offSet,
		bool polygonMesh, std::vector<Vertex>* barycenters) {
	MinimalAreaTriangulation<Real> MAT;
	std::vector<Point3D<Real> > vertices;
	std::vector<TriangleIndex> triangles;
	if(polygonMesh) {
		std::vector<CoredVertexIndex> vertices(edges.size());
		for(int i = 0; i != (int)edges.size(); ++i) {
			vertices[i].idx = edges[i].index;
			vertices[i].inCore = edges[i].inCore != 0;
		}
		mesh->addPolygon(vertices);
		return 1;
	}
	if(edges.size() > 3) {
		bool isCoplanar = false;

		if(barycenters) {
			for(unsigned i = 0; i != edges.size(); ++i) {
				for(unsigned j = 0; j != i; ++j) {
					if((i + 1) % edges.size() != j && (j + 1) % edges.size() != i) {
						Vertex v1 = edges[i].inCore ?
							mesh->inCorePoints(edges[i].index) :
							(*interiorVertices)[edges[i].index - offSet];
						Vertex v2 = edges[j].inCore ?
							mesh->inCorePoints(edges[j].index) :
							(*interiorVertices)[edges[j].index - offSet];
						for(int k = 0; k != 3; ++k) {
							if(v1.point[k] == v2.point[k]) isCoplanar = true;
						}
					}
				}
			}
		}
		if(isCoplanar) {
			Vertex c;
			c *= 0; // TODO: Ehm, what?
			for(int i = 0; i != (int)edges.size(); ++i) {
				c += edges[i].inCore ?
					mesh->inCorePoints(edges[i].index) :
					(*interiorVertices)[edges[i].index-offSet];
			}
			c /= (Real)edges.size();
			int cIdx = mesh->addOutOfCorePoint(c);
#pragma omp critical (add_barycenter_access)
			barycenters->push_back(c);
			for(int i = 0; i != (int)edges.size(); ++i) {
				std::vector<CoredVertexIndex> vertices(3);
				vertices[0].idx = edges[i].index;
				vertices[1].idx = edges[(i + 1) % edges.size()].index;
				vertices[2].idx = cIdx;
				vertices[0].inCore = edges[i].inCore != 0;
				vertices[1].inCore = edges[(i + 1) % edges.size()].inCore != 0;
				vertices[2].inCore = 0;
				mesh->addPolygon(vertices);
			}
			return edges.size();
		} else {
			vertices.resize(edges.size());
			// Add the points
			for(int i = 0; i != (int)edges.size(); ++i) {
				Vertex p = edges[i].inCore ?
					mesh->inCorePoints(edges[i].index) :
					(*interiorVertices)[edges[i].index-offSet];
				vertices[i] = p.point;
			}
			MAT.GetTriangulation(vertices, triangles);
			for(int i = 0; i != (int)triangles.size(); ++i) {
				std::vector<CoredVertexIndex> _vertices(3);
				for(int j = 0; j != 3; ++j) {
					_vertices[j].idx = edges[triangles[i].idx[j]].index;
					_vertices[j].inCore = edges[triangles[i].idx[j]].inCore != 0;
				}
				mesh->addPolygon(_vertices);
			}
		}
	} else if(edges.size() == 3) {
		std::vector<CoredVertexIndex> vertices(3);
		for(int i = 0; i != 3; ++i) {
			vertices[i].idx = edges[i].index;
			vertices[i].inCore = edges[i].inCore != 0;
		}
		mesh->addPolygon(vertices);
	}
	return (int)edges.size() - 2;
}

// TODO: Add voxel grid output test cases
template< int Degree , bool OutputDensity >
std::vector<Real> Octree<Degree, OutputDensity>::GetSolutionGrid(int& res, Real isoValue, int depth) {
	int maxDepth = boundaryType_ == BoundaryTypeNone ? tree_.maxDepth() - 1 : tree_.maxDepth();
	if(depth <= 0 || depth > maxDepth) depth = maxDepth;
	BSplineData<Degree, Real> fData;
	fData.set(boundaryType_ == BoundaryTypeNone ? depth + 1 : depth, boundaryType_);
	res = 1 << depth;
	fData.setValueTables();
	std::vector<Real> values(res * res * res);

	for(TreeOctNode* n = tree_.nextNode(); n; n = tree_.nextNode(n)) {
		if(n->depth() > (boundaryType_ == BoundaryTypeNone ? depth + 1 : depth)) continue;
		if(n->depth() < minDepth_) continue;
		int d;
		int idx[3];
		n->depthAndOffset(d, idx);
		int start[3];
		int end[3];
		for(int i = 0; i != 3; ++i) {
			// Get the index of the functions
			idx[i] = BinaryNode<double>::CenterIndex(d, idx[i]);
			// Figure out which samples fall into the range
			fData.setSampleSpan(idx[i], start[i], end[i]);
			// We only care about the odd indices
			if(!(start[i] & 1)) ++start[i];
			if(!(end[i] & 1)) --end[i];
			if(boundaryType_ == BoundaryTypeNone) {
				// (start[i]-1)>>1 >=   res/2 
				// (  end[i]-1)<<1 <  3*res/2
				start[i] = std::max(start[i], res + 1);
				end[i] = std::min(end[i], 3 * res - 1);
			}
		}
		Real coefficient = n->nodeData.solution;
		for(int x = start[0]; x <= end[0]; x += 2) {
			for(int y = start[1]; y <= end[1]; y += 2) {
				for(int z = start[2]; z <= end[2]; z += 2) {
					int xx = (x - 1) >> 1;
					int yy = (y - 1) >> 1;
					int zz = (z - 1) >> 1;
					if(boundaryType_ == BoundaryTypeNone) {
						xx -= res / 2;
						yy -= res / 2;
						zz -= res / 2;
					}
					values[zz * res * res + yy * res + xx] += coefficient *
						fData.valueTables(idx[0] + x * fData.functionCount()) *
						fData.valueTables(idx[1] + y * fData.functionCount()) *
						fData.valueTables(idx[2] + z * fData.functionCount());
				}
			}
		}
	}
	if(boundaryType_ == BoundaryTypeDirichlet)
		for(int i = 0; i != res * res * res; ++i) values[i] -= (Real)0.5;
	for(int i = 0; i != res * res * res; ++i) values[i] -= isoValue;

	return values;
}

////////////////
// VertexData //
////////////////
template<bool OutputDensity>
long long VertexData<OutputDensity>::CornerIndex(TreeOctNode const* node, int cIndex, int maxDepth) {
	int x[DIMENSION];
	int idx[DIMENSION];
	Cube::FactorCornerIndex(cIndex, x[0], x[1], x[2]);
	int d;
	int o[3];
	node->depthAndOffset(d, o);
	for(int i = 0; i != DIMENSION; ++i)
		idx[i] = BinaryNode<Real>::CornerIndex(maxDepth + 1, d, o[i], x[i]);
	return CornerIndexKey(idx);
}
