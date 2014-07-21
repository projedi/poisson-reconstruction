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

#pragma message( "[WARNING] Assuming that the number of octree nodes is less than INT_MAX" )

#ifndef NO_GRADIENT_DOMAIN_SOLUTION
#define GRADIENT_DOMAIN_SOLUTION 1
#else
#define GRADIENT_DOMAIN_SOLUTION 0
#endif
// Given the constraint vector-field V(p), there are two ways to solve for the coefficients, x, of the
// indicator function with respect to the B-spline basis {B_i(p)}
// 1] Find x minimizing:
//			|| V(p) - \sum_i \nabla x_i B_i(p) ||^2
//		which is solved by the system A_1x = b_1 where:
//			A_1[i,j] = < \nabla B_i(p) , \nabla B_j(p) >
//			b_1[i]   = < \nabla B_i(p) , V(p) >
// 2] Formulate this as a Poisson equation:
//			\sum_i x_i \Delta B_i(p) = \nabla \cdot V(p)
//		which is solved by the system A_2x = b_2 where:
//			A_2[i,j] = - < \Delta B_i(p) , B_j(p) >
//			b_2[i]   = - < B_i(p) , \nabla \cdot V(p) >
// Although the two system matrices should be the same (assuming that the B_i satisfy
// dirichlet/neumann boundary conditions) the constraint vectors can differ when V does not satisfy
// the Neumann boundary conditions:
//		A_1[i,j] = \int_R < \nabla B_i(p) , \nabla B_j(p) >
//               = \int_R [ \nabla \cdot ( B_i(p) \nabla B_j(p) ) - B_i(p) \Delta B_j(p) ]
//               = \int_dR < N(p) , B_i(p) \nabla B_j(p) > + A_2[i,j]
// and the first integral is zero if either f_i is zero on the boundary dR or the derivative of
// B_i across the boundary is zero. However, for the constraints we have:
//		b_1(i)   = \int_R < \nabla B_i(p) , V(p) >
//               = \int_R [ \nabla \cdot ( B_i(p) V(p) ) - B_i(p) \nabla \cdot V(p) ]
//               = \int_dR < N(p) ,  B_i(p) V(p) > + b_2[i]
// In particular, this implies that if the B_i satisfy the Neumann boundary conditions
// (rather than Dirichlet), and V is not zero across the boundary, then the two constraints are different.
// Forcing the < V(p) , N(p) > = 0 on the boundary, by killing off the component of the vector-field
// in the normal direction (FORCE_NEUMANN_FIELD), makes the two systems equal, and the value of this
// flag should be immaterial.
// Note that under interpretation 1, we have:
//		\sum_i b_1(i) = < \nabla \sum_ i B_i(p) , V(p) > = 0
// because the B_i's sum to one. However, in general, we could have
//		\sum_i b_2(i) \neq 0.
// This could cause trouble because the constant functions are in the kernel of the matrix A,
// so CG will misbehave if the constraint has a non-zero DC term.
// (Again, forcing < V(p) , N(p) > = 0 along the boundary resolves this problem.)

#ifndef NO_FORCE_NEUMANN_FIELD
#define FORCE_NEUMANN_FIELD 1
#else
#define FORCE_NEUMANN_FIELD 0
#endif
// This flag forces the normal component across the boundary of the integration domain to be zero.
// This should be enabled if GRADIENT_DOMAIN_SOLUTION is not, so that CG doesn't run into trouble.

// 0 produces segmentation fault on the original code
#ifndef SPLAT_ORDER_1
#define SPLAT_ORDER 2
#else
#define SPLAT_ORDER 1
#endif

size_t const MEMORY_ALLOCATOR_BLOCK_SIZE = 1 << 12;

#if !FORCE_NEUMANN_FIELD
#pragma message("[WARNING] Not zeroing out normal component on boundary")
#endif

#include "BSplineData.h"
#include "HashMap.h"
#include "Octree.h"
#include "PPolynomial.h"
#include "Ply.h"
#include "SparseMatrix.h"
#include "Time.h"

typedef float Real;
typedef float MatrixReal;

template<bool StoreDensity>
class TreeNodeData {
public:
	int nodeIndex;
	union {
		int mcIndex;
		int normalIndex;
	};
	Real centerWeightContribution[StoreDensity ? 2 : 1];
	Real constraint;
	Real solution;
	int pointIndex;

	TreeNodeData();
};

template<bool OutputDensity>
class RootInfo {
	typedef OctNode<TreeNodeData<OutputDensity>, Real> TreeOctNode;
public:
	TreeOctNode const* node;
	int edgeIndex;
	long long key;
};

template<bool OutputDensity>
class VertexData {
public:
	typedef OctNode<TreeNodeData<OutputDensity>, Real> TreeOctNode;
	static long long CornerIndex(TreeOctNode const* node, int cIndex, int maxDepth);
private:
	static int const VERTEX_COORDINATE_SHIFT = (sizeof(long long) * 8) / 3;
	static long long CornerIndexKey(int const idx[DIMENSION])
		{ return (long long)idx[0] | ((long long)idx[1]) << VERTEX_COORDINATE_SHIFT |
			((long long)idx[2]) << (2 * VERTEX_COORDINATE_SHIFT); }
};

template<int Degree>
struct Indices {
	Indices() { memset(idx, -1, sizeof(int) * Degree); }
	int& operator[](int i) { return idx[i]; }
	int operator[](int i) const { return idx[i]; }
private:
	int idx[Degree];
};

template<class TreeOctNode, class IndicesS>
struct TableData {
	TableData(): count_(0) { }
	virtual ~TableData() { clear(); }
	void clear() { table_.clear(); count_ = 0; }

	IndicesS& operator[](TreeOctNode const* node) { return indices(node); }
	IndicesS const& operator[](TreeOctNode const* node) const { return indices(node); }

	IndicesS& indices(TreeOctNode const* node)
		{ return table_[node->nodeData.nodeIndex + offsets_[node->depth()]]; }
	IndicesS const& indices(TreeOctNode const* node) const
		{ return table_[node->nodeData.nodeIndex + offsets_[node->depth()]]; }

	int count() const { return count_; }

	int& offsets(int i) { return offsets_[i]; }
	int offsets(int i) const { return offsets_[i]; }

	void setCount(int count) { count_ = count; }
	void resizeOffsets(int size, int val) { offsets_.resize(size, val); }
	void resizeTable(int size) { table_.resize(size); }
private:
	int count_;
	std::vector<IndicesS> table_;
	std::vector<int> offsets_;
};

template<bool OutputDensity>
class SortedTreeNodes {
public:
	typedef OctNode<TreeNodeData<OutputDensity>, Real> TreeOctNode;
	typedef typename TreeOctNode::ConstNeighborKey3 TreeConstNeighborKey3;
	typedef Indices<Cube::CORNERS> CornerIndices;
	typedef TableData<TreeOctNode, CornerIndices> CornerTableData;
	typedef Indices<Cube::EDGES> EdgeIndices;
	typedef TableData<TreeOctNode, EdgeIndices> EdgeTableData;

	std::vector<int> nodeCount;
	std::vector<TreeOctNode*> treeNodes;
	int maxDepth;

	SortedTreeNodes(): maxDepth(0) { }
	void set(TreeOctNode& root);

// TODO: setTable and getMaxCount between Corner and Edge share a lot of code. But straight up
// TODO: extraction only makes it worse. Refactor it somehow.

	void setCornerTable(CornerTableData& cData, TreeOctNode const* rootNode, int depth, int threads) const;
	void setCornerTable(CornerTableData& cData, TreeOctNode const* rootNode, int threads) const
		{ setCornerTable(cData, rootNode, maxDepth - 1, threads); }
	int getMaxCornerCount(int depth, int maxDepth, int threads) const;

	void setEdgeTable(EdgeTableData& eData, TreeOctNode const* rootNode, int depth, int threads);
	void setEdgeTable(EdgeTableData& eData, TreeOctNode const* rootNode, int threads)
		{ setEdgeTable(eData, rootNode, maxDepth - 1, threads); }
	int getMaxEdgeCount(TreeOctNode const* rootNode, int depth, int threads) const;
};

struct PointData {
	Point3D<Real> position;
	Real weight;
	Real coarserValue;
	PointData(Point3D<Real> p = Point3D<Real>(), Real w = 0): position(p), weight(w), coarserValue(0) { }
};

template<class C, int N>
class Stencil {
public:
	C& at(int i, int j, int k) { return values[i][j][k]; }
	C const& at(int i, int j, int k) const { return values[i][j][k]; }
private:
	C values[N][N][N];
};

typedef Stencil<Point3D<double>, 5> DivergenceStencil;
typedef Stencil<DivergenceStencil, 2> DivergenceStencils;

typedef Stencil<double, 5> LaplacianStencil;
typedef Stencil<LaplacianStencil, 2> LaplacianStencils;

typedef Stencil<double, 3> CenterEvaluationStencil;
typedef Stencil<CenterEvaluationStencil, 2> CenterEvaluationStencils;

typedef Stencil<Stencil<double, 3>, 2> CornerEvaluationStencil;
typedef Stencil<CornerEvaluationStencil, 2> CornerEvaluationStencils;

typedef Stencil<Stencil<Point3D<double>, 5>, 2> CornerNormalEvaluationStencil;
typedef Stencil<CornerNormalEvaluationStencil, 2> CornerNormalEvaluationStencils;

struct CenterValueStencil {
	CenterEvaluationStencil stencil;
	CenterEvaluationStencils stencils;
};

struct CornerValueStencil {
	CornerEvaluationStencil stencil;
	CornerEvaluationStencils stencils;
};

struct CornerNormalStencil {
	CornerNormalEvaluationStencil stencil;
	CornerNormalEvaluationStencils stencils;
};

struct UpSampleData {
	UpSampleData(): start(0) {
		v[0] = 0;
		v[1] = 0;
	}
	UpSampleData(int s, double v1, double v2): start(s) {
		v[0] = v1;
		v[1] = v2;
	}
	int start;
	double v[2];
};

// For computing the iso-surface there is a lot of re-computation of information across shared geometry.
// For function values we don't care so much.
// For edges we need to be careful so that the mesh remains water-tight
template<bool OutputDensity>
struct RootData: SortedTreeNodes<OutputDensity>::CornerTableData,
		SortedTreeNodes<OutputDensity>::EdgeTableData {
	// Edge to iso-vertex map
	HashMap<long long, int> boundaryRoots;
	// Vertex to ( value , normal ) map
	HashMap<long long, std::pair<Real, Point3D<Real> > > boundaryValues;

	std::vector<int> interiorRoots;
	std::vector<Real> cornerValues;
	std::vector<Point3D<Real> > cornerNormals;
	std::vector<char> cornerValuesSet;
	std::vector<char> cornerNormalsSet;
	std::vector<char> edgesSet;

	int cCount() const { return SortedTreeNodes<OutputDensity>::CornerTableData::count(); }
	int eCount() const { return SortedTreeNodes<OutputDensity>::EdgeTableData::count(); }

	int cornerIndices(typename SortedTreeNodes<OutputDensity>::TreeOctNode const* node, int idx)
		{ return SortedTreeNodes<OutputDensity>::CornerTableData::indices(node)[idx]; }
	int edgeIndices(typename SortedTreeNodes<OutputDensity>::TreeOctNode const* node, int idx)
		{ return SortedTreeNodes<OutputDensity>::EdgeTableData::indices(node)[idx]; }
};

struct Range3D {
	static Range3D FullRange() {
		Range3D range;
		range.xStart = range.yStart = range.zStart = 0;
		range.xEnd = range.yEnd = range.zEnd = 5;
		return range;
	}

	int xStart;
	int xEnd;
	int yStart;
	int yEnd;
	int zStart;
	int zEnd;
};

template<int Degree, bool OutputDensity>
class Octree {
public:
	typedef OctNode<TreeNodeData<OutputDensity>, Real> TreeOctNode;
	typedef typename TreeOctNode::NeighborKey3 TreeNeighborKey3;
	typedef typename TreeOctNode::ConstNeighborKey3 TreeConstNeighborKey3;

	static double maxMemoryUsage() { return ((double)maxMemoryUsage_) / (1 << 20); }
	static void resetMaxMemoryUsage() { maxMemoryUsage_ = 0; }

	Octree(int threads, int maxDepth, BoundaryType boundaryType);

	void finalize(int subdivisionDepth);
	std::vector<Real> GetSolutionGrid(int& res, Real isoValue, int depth);
	int setTree(std::string const& fileName, int maxDepth, int minDepth, int kernelDepth, Real samplesPerNode,
		Real scaleFactor, bool useConfidence, bool useNormalWeights, Real constraintWeight,
		int adaptiveExponent, XForm<Real, 4> xForm);

	void SetLaplacianConstraints();
	void ClipTree();
	int LaplacianMatrixIteration(int subdivideDepth, bool showResidual, int minIters, double accuracy,
			int maxSolveDepth, int fixedIters);

	Real GetIsoValue() const;
	template<class Vertex>
	void GetMCIsoTriangles(Real isoValue, int subdivideDepth, CoredFileMeshData<Vertex>* mesh,
			int nonLinearFit, bool addBarycenter, bool polygonMesh);

	TreeOctNode const& tree() const { return tree_; }
private:
	typedef typename BSplineData<Degree, Real>::Integrator Integrator;
	typedef typename TreeOctNode::Neighbors3 TreeNeighbors3;
	typedef typename TreeOctNode::Neighbors5 TreeNeighbors5;
	typedef typename TreeOctNode::ConstNeighbors3 TreeConstNeighbors3;
	typedef typename TreeOctNode::ConstNeighbors5 TreeConstNeighbors5;
	typedef typename BSplineData<Degree, Real>::template CenterEvaluator<1> CenterEvaluator1;
	typedef typename BSplineData<Degree, Real>::template CornerEvaluator<2> CornerEvaluator2;
	typedef typename SortedTreeNodes<OutputDensity>::CornerTableData CornerTableData;
	typedef std::pair<RootInfo<OutputDensity>, RootInfo<OutputDensity> > edge_t;
	typedef std::vector<edge_t> edges_t;
	typedef HashMap<long long, std::pair<RootInfo<OutputDensity>, int> > vertex_count_t;

	class FaceEdgesFunction {
	public:
		FaceEdgesFunction(int maxDepth, edges_t& edges, vertex_count_t& vertexCount,
				TreeConstNeighborKey3& neighborKey3):
			maxDepth(maxDepth), edges(edges), vertexCount(vertexCount), neighborKey3(neighborKey3) { }
		void setFIndex(int fIndex) { this->fIndex = fIndex; }
		void operator()(TreeOctNode const* node1, TreeOctNode const* node2) const;
	private:
		int fIndex;
		int maxDepth;
		edges_t& edges;
		vertex_count_t& vertexCount;
		TreeConstNeighborKey3& neighborKey3;
	};

	class SetDivergenceStencilFunction {
	public:
		SetDivergenceStencilFunction(Octree const& o, int depth, Integrator const& integrator,
				bool scatter, int center):
			o(o), depth(depth), integrator(integrator), scatter(scatter), center(center) { }
		Point3D<double> operator()(int x, int y, int z) const {
			int offset[] = { center, center, center };
			int _offset[] = { x + center - 2, y + center - 2, z + center - 2 };
			return scatter ?
				o.GetDivergence1(integrator, depth, offset, _offset, false) :
				o.GetDivergence2(integrator, depth, offset, _offset, false);
		}
	private:
		Octree const& o;
		int depth;
		Integrator const& integrator;
		bool scatter;
		int center;
	};

	class SetDivergenceStencilsFunction {
	public:
		SetDivergenceStencilsFunction(Octree const& o, int depth, Integrator const& integrator,
				bool scatter, int center):
			o(o), depth(depth), integrator(integrator), scatter(scatter), center(center) { }
		Point3D<double> operator()(int i, int j, int k, int x, int y, int z) const {
			int offset[] = { center + i, center + j, center + k };
			int _offset[] = { x + center / 2 - 2, y + center / 2 - 2, z + center / 2 - 2 };
			return scatter ?
				o.GetDivergence1(integrator, depth, offset, _offset, true) :
				o.GetDivergence2(integrator, depth, offset, _offset, true);
		}
	private:
		Octree const& o;
		int depth;
		Integrator const& integrator;
		bool scatter;
		int center;
	};

	class SetLaplacianStencilFunction {
	public:
		SetLaplacianStencilFunction(Octree const& o, int depth, Integrator const& integrator, int center):
			o(o), depth(depth), integrator(integrator), center(center) { }
		double operator()(int x, int y, int z) const {
			int offset[] = { center, center, center };
			int _offset[] = { x + center - 2, y + center - 2, z + center - 2 };
			return o.GetLaplacian(integrator, depth, offset, _offset, false);
		}
	private:
		Octree const& o;
		int depth;
		Integrator const& integrator;
		int center;
	};

	class SetLaplacianStencilsFunction {
	public:
		SetLaplacianStencilsFunction(Octree const& o, int depth, Integrator const& integrator, int center):
			o(o), depth(depth), integrator(integrator), center(center) { }
		double operator()(int i, int j, int k, int x, int y, int z) const {
			int offset[] = { center + i, center + j, center + k };
			int _offset[] = { x + center / 2 - 2, y + center / 2 - 2, z + center / 2 - 2 };
			return o.GetLaplacian(integrator, depth, offset, _offset, true);
		}
	private:
		Octree const& o;
		int depth;
		Integrator const& integrator;
		int center;
	};

	class SetCenterEvaluationStencilFunction {
	public:
		SetCenterEvaluationStencilFunction(int depth, int center, CenterEvaluator1 const& evaluator):
			depth(depth), center(center), evaluator(evaluator) { }
		double operator()(int x, int y, int z) const {
			int _offset[] = { x + center - 1, y + center - 1, z + center - 1 };
			return evaluator.value(depth, center, _offset[0], false, false) *
				evaluator.value(depth, center, _offset[1], false, false) *
				evaluator.value(depth, center, _offset[2], false, false);
		}
	private:
		int depth;
		int center;
		CenterEvaluator1 const& evaluator;
	};

	class SetCenterEvaluationStencilsFunction {
	public:
		SetCenterEvaluationStencilsFunction(int depth, int center, CenterEvaluator1 const& evaluator):
			depth(depth), center(center), evaluator(evaluator) { }
		double operator()(int cx, int cy, int cz, int x, int y, int z) const {
			int idx[] = { center + cx, center + cy, center + cz };
			int off[] = { x + center / 2 - 1, y + center / 2 - 1, z + center / 2 - 1};
			return evaluator.value(depth, idx[0], off[0], false, true) *
				evaluator.value(depth, idx[1], off[1], false, true) *
				evaluator.value(depth, idx[2], off[2], false, true);
		}
	private:
		int depth;
		int center;
		CenterEvaluator1 const& evaluator;
	};

	class SetCornerEvaluationStencilFunction {
	public:
		SetCornerEvaluationStencilFunction(int depth, int center, CornerEvaluator2 const& evaluator):
			depth(depth), center(center), evaluator(evaluator) { }
		double operator()(int cx, int cy, int cz, int x, int y, int z) const {
			int off[] = { center + x - 1, center + y - 1, center + z - 1 };
			return evaluator.value(depth, center, cx, off[0], false, false) *
				evaluator.value(depth, center, cy, off[1], false, false) *
				evaluator.value(depth, center, cz, off[2], false, false);
		}
	private:
		int depth;
		int center;
		CornerEvaluator2 const& evaluator;
	};

	class SetCornerEvaluationStencilsFunction {
	public:
		SetCornerEvaluationStencilsFunction(int depth, int center, CornerEvaluator2 const& evaluator):
			depth(depth), center(center), evaluator(evaluator) { }
		double operator()(int cx, int cy, int cz, int _cx, int _cy, int _cz,
				int x, int y, int z) const {
			int idx[] = { center + _cx, center + _cy, center + _cz };
			int off[] = { center / 2 + x - 1, center / 2 + y - 1, center / 2 + z - 1 };
			return evaluator.value(depth, idx[0], cx, off[0], false, true) *
				evaluator.value(depth, idx[1], cy, off[1], false, true) *
				evaluator.value(depth, idx[2], cz, off[2], false, true);
		}
	private:
		int depth;
		int center;
		CornerEvaluator2 const& evaluator;
	};

	class SetCornerNormalEvaluationStencilFunction {
	public:
		SetCornerNormalEvaluationStencilFunction(int depth, int center, CornerEvaluator2 const& evaluator):
			depth(depth), center(center), evaluator(evaluator) { }
		Point3D<double> operator()(int cx, int cy, int cz, int x, int y, int z) const {
			int off[] = { center + x - 2, center + y - 2, center + z - 2 };
			double v[] = { evaluator.value(depth, center, cx, off[0], false, false),
				evaluator.value(depth, center, cy, off[1], false, false),
				evaluator.value(depth, center, cz, off[2], false, false) };
			double dv[] = { evaluator.value(depth, center, cx, off[0], true, false),
				evaluator.value(depth, center, cy, off[1], true, false),
				evaluator.value(depth, center, cz, off[2], true, false) };
			return Point3D<double>(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]);
		}
	private:
		int depth;
		int center;
		CornerEvaluator2 const& evaluator;
	};

	class SetCornerNormalEvaluationStencilsFunction {
	public:
		SetCornerNormalEvaluationStencilsFunction(int depth, int center, CornerEvaluator2 const& evaluator):
			depth(depth), center(center), evaluator(evaluator) { }
		Point3D<double> operator()(int cx, int cy, int cz, int _cx, int _cy, int _cz,
				int x, int y, int z) const {
			int idx[] = { center + _cx, center + _cy, center + _cz };
			int off[] = { center / 2 + x - 2, center / 2 + y - 2, center / 2 + z - 2 };
			double v[] = { evaluator.value(depth, idx[0], cx, off[0], false, true),
				evaluator.value(depth, idx[1], cy, off[1], false, true),
				evaluator.value(depth, idx[2], cz, off[2], false, true) };
			double dv[] = { evaluator.value(depth, idx[0], cx, off[0], true, true),
				evaluator.value(depth, idx[1], cy, off[1], true, true),
				evaluator.value(depth, idx[2], cz, off[2], true, true) };
			return Point3D<double>(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]);
		}
	private:
		int depth;
		int center;
		CornerEvaluator2 const& evaluator;
	};

	class UpSampleCoarserSolutionFunction {
	public:
		UpSampleCoarserSolutionFunction(Vector<Real>& Solution, size_t start):
			Solution(Solution), start(start) { }
		void operator()(int i, TreeOctNode const* node, UpSampleData* usData, int* idxs) const {
			double dxyz = usData[0].v[idxs[0]] * usData[1].v[idxs[1]] * usData[2].v[idxs[2]];
			Solution[i - start] += (Real)(node->nodeData.solution * dxyz);
		}
	private:
		Vector<Real>& Solution;
		size_t start;
	};

	template<class C>
	class DownSampleFunction {
	public:
		DownSampleFunction(C* constraints): constraints(constraints) { }
		void operator()(int i, TreeOctNode const* node, UpSampleData* usData, int* idxs) const {
			C cx = constraints[i] * usData[0].v[idxs[0]];
			C cxy = cx * usData[1].v[idxs[1]];
			C cxyz = cxy * usData[2].v[idxs[2]];
#pragma omp atomic
			constraints[node->nodeData.nodeIndex] += cxyz;
		}
	private:
		C* constraints;
	};

	template<class C>
	class UpSample1Function {
	public:
		UpSample1Function(C* coefficients): coefficients(coefficients) { }
		void operator()(int i, TreeOctNode const* node, UpSampleData* usData, int* idxs) const {
			double dx = usData[0].v[idxs[0]];
			double dxy = dx * usData[1].v[idxs[1]];
			double dxyz = dxy * usData[2].v[idxs[2]];
			coefficients[i] += coefficients[node->nodeData.nodeIndex] * (Real)dxyz;
		}
	private:
		C* coefficients;
	};

	template<class C>
	class UpSample2Function {
	public:
		UpSample2Function(int depth, C* fineCoefficients, C const* coarseCoefficients,
				SortedTreeNodes<OutputDensity> const& sNodes): depth(depth),
			fineCoefficients(fineCoefficients), coarseCoefficients(coarseCoefficients), sNodes(sNodes) { }
		void operator()(int i, TreeOctNode const* node, UpSampleData* usData, int* idxs) const {
			double dx = usData[0].v[idxs[0]];
			double dxy = dx * usData[1].v[idxs[1]];
			double dxyz = dxy * usData[2].v[idxs[2]];
			fineCoefficients[i - sNodes.nodeCount[depth]] +=
				coarseCoefficients[node->nodeData.nodeIndex - sNodes.nodeCount[depth - 1]] * (Real)dxyz;
		}
	private:
		int depth;
		C* fineCoefficients;
		C const* coarseCoefficients;
		SortedTreeNodes<OutputDensity> const& sNodes;
	};

	class GetFixedDepthLaplacianGetNodeFunction {
	public:
		GetFixedDepthLaplacianGetNodeFunction(SortedTreeNodes<OutputDensity> const& sNodes, size_t start):
			sNodes(sNodes), start(start) { }
		TreeOctNode* operator()(int i) const {
			return sNodes.treeNodes[i + start];
		}
	private:
		SortedTreeNodes<OutputDensity> const& sNodes;
		size_t start;
	};

	class GetFixedDepthLaplacianGetRowSizeFunction {
	public:
		GetFixedDepthLaplacianGetRowSizeFunction(Octree& o): o(o) { }
		int operator()(TreeNeighbors5 const& neighbors5, bool symmetric) const {
			return o.GetMatrixRowSize(neighbors5, symmetric);
		}
	private:
		Octree& o;
	};

	class GetFixedDepthLaplacianSetRowFunction {
	public:
		GetFixedDepthLaplacianSetRowFunction(Octree& o): o(o) { }
		int operator()(TreeNeighbors5 const& neighbors5, SparseSymmetricMatrix<MatrixReal>& m, int row,
				int offset, Integrator const& integrator, Stencil<double, 5> const& stencil,
				bool symmetric) const {
			return o.SetMatrixRow(neighbors5, m, row, offset, integrator, stencil, symmetric);
		}
	private:
		Octree& o;
	};

	class GetRestrictedFixedDepthLaplacianGetNodeFunction {
	public:
		GetRestrictedFixedDepthLaplacianGetNodeFunction(Octree& o,
				SortedTreeNodes<OutputDensity> const& sNodes, int depth, std::vector<int> const& entries,
				int rDepth, int rOff[3], std::vector<Range3D>& ranges):
			o(o), sNodes(sNodes), depth(depth), entries(entries), rDepth(rDepth), rOff(rOff), ranges(ranges) { }
		TreeOctNode* operator()(int i) const {
			TreeOctNode* node = sNodes.treeNodes[entries[i]];
			int d;
			int off[3];
			node->depthAndOffset(d, off);
			off[0] >>= depth - rDepth;
			off[1] >>= depth - rDepth;
			off[2] >>= depth - rDepth;
			bool isInterior = off[0] == rOff[0] && off[1] == rOff[1] && off[2] == rOff[2];

			if(!isInterior) o.SetMatrixRowBounds(node, rDepth, rOff, ranges[omp_get_thread_num()]);
			else ranges[omp_get_thread_num()] = Range3D::FullRange();
			return node;
		}
	private:
		Octree& o;
		SortedTreeNodes<OutputDensity> const& sNodes;
		int depth;
		std::vector<int> const& entries;
		int rDepth;
		int* rOff;
		std::vector<Range3D>& ranges;
	};

	class GetRestrictedFixedDepthLaplacianGetRowSizeFunction {
	public:
		GetRestrictedFixedDepthLaplacianGetRowSizeFunction(Octree& o, std::vector<Range3D> const& ranges):
			o(o), ranges(ranges) { }
		int operator()(TreeNeighbors5 const& neighbors5, bool symmetric) const {
			return o.GetMatrixRowSize(neighbors5, ranges[omp_get_thread_num()], symmetric);
		}
	private:
		Octree& o;
		std::vector<Range3D> const& ranges;
	};

	class GetRestrictedFixedDepthLaplacianSetRowFunction {
	public:
		GetRestrictedFixedDepthLaplacianSetRowFunction(Octree& o, std::vector<Range3D> const& ranges):
			o(o), ranges(ranges) { }
		int operator()(TreeNeighbors5 const& neighbors5, SparseSymmetricMatrix<MatrixReal>& m, int row, int,
				Integrator const& integrator, Stencil<double, 5> const& stencil, bool symmetric) const {
			return o.SetMatrixRow(neighbors5, m, row, 0, integrator, stencil,
					ranges[omp_get_thread_num()], symmetric);
		}
	private:
		Octree& o;
		std::vector<Range3D> const& ranges;
	};

	class SplatOrientedPointGetNeighborsFunction {
	public:
		SplatOrientedPointGetNeighborsFunction(TreeNeighborKey3& key): neighborKey(key) { }
		TreeConstNeighbors3& operator()(TreeOctNode* node) const {
			return (TreeConstNeighbors3&)neighborKey.setNeighbors(node);
		}
	private:
		TreeNeighborKey3& neighborKey;
	};

	class GetRootGetNeighborsFunction {
	public:
		GetRootGetNeighborsFunction(TreeConstNeighborKey3& key): neighborKey3(key) { }
		TreeConstNeighbors3& operator()(TreeOctNode const* node) const {
			return neighborKey3.getNeighbors3(node);
		}
	private:
		TreeConstNeighborKey3& neighborKey3;
	};

	static double MemoryUsage();
	static void UpdateCoarserSupportBounds(TreeOctNode const* node, Range3D& range);
	static int IsBoundaryFace(TreeOctNode const* node, int faceIndex, int subdivideDepth);
	static int IsBoundaryEdge(TreeOctNode const* node, int edgeIndex, int subdivideDepth);
	static int IsBoundaryEdge(TreeOctNode const* node, int dir, int x, int y, int subidivideDepth);
	template<class Vertex>
	static int AddTriangles(CoredFileMeshData<Vertex>* mesh, std::vector<CoredPointIndex>& edges,
			std::vector<Vertex>* interiorVertices, int offSet, bool polygonMesh,
			std::vector<Vertex>* barycenters);
	static std::vector<edges_t> GetEdgeLoops(edges_t& edges);
	static int GetRootIndex(TreeOctNode const* node, int edgeIndex, int maxDepth,
			TreeConstNeighborKey3& neighborKey3, RootInfo<OutputDensity>& ri);
	static int GetRootIndex(RootInfo<OutputDensity> const& ri, RootData<OutputDensity>& rootData,
			CoredPointIndex& index);
	static int GetRootPair(RootInfo<OutputDensity> const& root, int maxDepth,
			TreeConstNeighborKey3& neighborKey3, RootInfo<OutputDensity>& pair);
	static bool IsInset(TreeOctNode const* node);

	int refineBoundary(int subdivisionDepth);
	bool inBounds(Point3D<Real>) const;
	double GetLaplacian(Integrator const& integrator, int d, int const off1[3], int const off2[3],
			bool childParent) const;
	double GetDivergence1(Integrator const& integrator, int d, int const off1[3], int const off2[3],
			bool childParent, Point3D<Real> const& normal1) const;
	Point3D<double> GetDivergence1(Integrator const& integrator, int d, int const off1[3],
			int const off2[3], bool childParent) const;
	double GetDivergence2(Integrator const& integrator, int d, int const off1[3], int const off2[3],
			bool childParent, Point3D<Real> const& normal2) const;
	Point3D<double> GetDivergence2(Integrator const& integrator, int d, int const off1[3],
			int const off2[3], bool childParent) const;
	int SolveFixedDepthMatrix(int depth, Integrator const& integrator,
			SortedTreeNodes<OutputDensity> const& sNodes, Real* subConstraints,
			bool showResidual, int minIters, double accuracy, bool noSolve, int fixedIters);
	int SolveFixedDepthMatrix(int depth, Integrator const& integrator,
			SortedTreeNodes<OutputDensity> const& sNodes, Real* subConstraints, int startingDepth,
			bool showResidual, int minIters, double accuracy, bool noSolve, int fixedIters);
	void SetMatrixRowBounds(TreeOctNode const* node, int rDepth, int const rOff[3], 
			Range3D& range) const;
	int GetMatrixRowSize(TreeNeighbors5 const& neighbors5, bool symmetric) const
		{ return GetMatrixRowSize(neighbors5, Range3D::FullRange(), symmetric); }
	int GetMatrixRowSize(TreeNeighbors5 const& neighbors5, Range3D const& range, bool symmetric) const;
	int SetMatrixRow(TreeNeighbors5 const& neighbors5, SparseSymmetricMatrix<MatrixReal>& m, int row,
			int off, Integrator const& integrator, Stencil<double, 5> const& stencil,
			bool symmetric) const
		{ return SetMatrixRow(neighbors5, m, row, off, integrator, stencil, Range3D::FullRange(), symmetric); }
	int SetMatrixRow(TreeNeighbors5 const& neighbors5, SparseSymmetricMatrix<MatrixReal>& m, int row,
			int offset, Integrator const& integrator, Stencil<double, 5> const& stencil,
			Range3D const& range, bool symmetric) const;
	LaplacianStencil SetLaplacianStencil(int depth, Integrator const& integrator) const;
	LaplacianStencils SetLaplacianStencils(int depth, Integrator const& integrator) const;
	DivergenceStencil SetDivergenceStencil(int depth, Integrator const& integrator, bool scatter) const;
	DivergenceStencils SetDivergenceStencils(int depth, Integrator const& integrator, bool scatter) const;
	CenterEvaluationStencil SetCenterEvaluationStencil(CenterEvaluator1 const& evaluator, int depth) const;
	CenterEvaluationStencils SetCenterEvaluationStencils(CenterEvaluator1 const& evaluator, int depth) const;
	CornerEvaluationStencil SetCornerEvaluationStencil(CornerEvaluator2 const& evaluator, int depth) const;
	CornerEvaluationStencils SetCornerEvaluationStencils(CornerEvaluator2 const& evaluator, int depth) const;
	CornerNormalEvaluationStencil SetCornerNormalEvaluationStencil(CornerEvaluator2 const& evaluator,
			int depth) const;
	CornerNormalEvaluationStencils SetCornerNormalEvaluationStencils(CornerEvaluator2 const& evaluator,
			int depth) const;
	void UpdateConstraintsFromCoarser(TreeNeighbors5 const& neighbors5, TreeNeighbors5 const& pNeighbors5,
			TreeOctNode* node, Real const* metSolution, Integrator const& integrator,
			Stencil<double, 5> const& stencil) const;
	void SetCoarserPointValues(int depth, SortedTreeNodes<OutputDensity> const& sNodes, Real* metSolution);
	Real WeightedCoarserFunctionValue(TreeNeighborKey3 const& neighborKey3, TreeOctNode const* node,
			Real* metSolution) const;
	Vector<Real> UpSampleCoarserSolution(int depth, SortedTreeNodes<OutputDensity> const& sNodes) const;
	template<class C>
	void DownSample(int depth, SortedTreeNodes<OutputDensity> const& sNodes, C* constraints) const;
	template<class C>
	void UpSample(int depth, SortedTreeNodes<OutputDensity> const& sNodes, C* coefficients) const;
	template<class C>
	void UpSample(int depth, SortedTreeNodes<OutputDensity> const& sNodes, C const* coarseCoefficients,
			C* fineCoefficients) const;
	template<class F1, class F2, class F3>
	SparseSymmetricMatrix<Real> GetFixedDepthLaplacianGeneric(int depth, Integrator const& integrator,
			SortedTreeNodes<OutputDensity> const& sNodes, Real const* metSolution, size_t range,
			F1 const& getNode, F2 const& getRowSize, F3 const& setRow);
	SparseSymmetricMatrix<Real> GetFixedDepthLaplacian(int depth, Integrator const& integrator,
			SortedTreeNodes<OutputDensity> const& sNodes, Real const* metSolution);
	SparseSymmetricMatrix<Real> GetRestrictedFixedDepthLaplacian(int depth, Integrator const& integrator,
			std::vector<int> const& entries, int entryCount, TreeOctNode const* rNode, Real radius,
			SortedTreeNodes<OutputDensity> const& sNodes, Real const* metSolution);
	void SetIsoCorners(Real isoValue, TreeOctNode* leaf, CornerTableData& cData, char* valuesSet,
			Real* values, TreeConstNeighborKey3& nKey, std::vector<Real> const& metSolution,
			CornerEvaluator2 const& evaluator, CornerEvaluationStencil const&, CornerEvaluationStencils const&);
	template<class Vertex>
	int SetMCRootPositions(TreeOctNode* node, int sDepth, Real isoValue,
			TreeConstNeighborKey3& neighborKey3, RootData<OutputDensity>& rootData,
			std::vector<Vertex>* interiorVertices, CoredFileMeshData<Vertex>* mesh,
			std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
			CornerNormalEvaluationStencil const&, CornerNormalEvaluationStencils const&, bool nonLinearFit);
	template<class Vertex>
	int GetMCIsoTriangles(TreeOctNode* node, TreeConstNeighborKey3& neighborKey3,
			CoredFileMeshData<Vertex>* mesh, RootData<OutputDensity>& rootData,
			std::vector<Vertex>* interiorVertices, int offSet, int sDepth, bool polygonMesh,
			std::vector<Vertex>* barycenters);
	void GetMCIsoEdges(TreeOctNode* node, TreeConstNeighborKey3& neighborKey3, int sDepth, edges_t& edges);
	template<class Vertex>
	int GetRoot(RootInfo<OutputDensity> const& ri, Real isoValue, TreeConstNeighborKey3& neighborKey3,
			Vertex& vertex, RootData<OutputDensity>& rootData, int sDepth, std::vector<Real> const& metSolution,
			CornerEvaluator2 const& evaluator, CornerNormalEvaluationStencil const&,
			CornerNormalEvaluationStencils const&, bool nonLinearFit);
	void UpdateWeightContribution(TreeOctNode* node, Point3D<Real> const& position,
			TreeNeighborKey3& neighborKey, Real weight = 1.0) const;
	Real GetSampleWeight(TreeOctNode const* node, Point3D<Real> const& position,
			TreeConstNeighbors3& neighbors) const;
	// To abstract from "TreeOctNode" and "TreeOctNode const"
	template<class OctNodeS, class GetNeighbors>
	void GetSampleDepthAndWeight(OctNodeS* node, Point3D<Real> const& position,
			GetNeighbors const& getNeighbors, Real samplesPerNode, Real& depth, Real& weight) const;
	void SplatOrientedPoint(TreeOctNode* node, Point3D<Real> const& point, Point3D<Real> const& normal,
			TreeNeighborKey3& neighborKey);
	Real SplatOrientedPoint(Point3D<Real> const& point, Point3D<Real> const& normal,
			TreeNeighborKey3& neighborKey, int kernelDepth, Real samplesPerNode, int minDepth, int maxDepth);

	bool HasNormals(TreeOctNode* node, Real epsilon) const;
	Real getCenterValue(TreeConstNeighborKey3 const& neighborKey3, TreeOctNode const* node,
			std::vector<Real> const& metSolution, CenterEvaluator1 const& evaluator,
			Stencil<double, 3> const& stencil, Stencil<double, 3> const& pStencil, bool isInterior) const;
	Real getCornerValue(TreeConstNeighborKey3 const& neighborKey3, TreeOctNode const* node, int corner,
			std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
			Stencil<double, 3> const& stencil, CornerEvaluationStencil const& stencils, bool isInterior) const;
	Point3D<Real> getCornerNormal(TreeConstNeighbors5 const& neighbors5,
			TreeConstNeighbors5 const& pNeighbors5, TreeOctNode const* node, int corner,
			std::vector<Real> const& metSolution, CornerEvaluator2 const& evaluator,
			Stencil<Point3D<double>, 5> const& nStencil, CornerNormalEvaluationStencil const&,
			bool isInterior) const;
private:
	static size_t maxMemoryUsage_;

	int threads_;
	BoundaryType boundaryType_;
	Real radius_;
	int width_;
	Real postDerivativeSmooth_;
	bool constrainValues_;
	TreeOctNode tree_;
	std::vector<Point3D<Real> > normals_;
	BSplineData<Degree, Real> fData_;
	SortedTreeNodes<OutputDensity> sNodes_;
	Real samplesPerNode_;
	int splatDepth_;
	int minDepth_;
	Real scale_;
	Point3D<Real> center_;
	std::vector<PointData> points_;
};

#include "MultiGridOctreeData.inl"
