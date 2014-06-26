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

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <omp.h>

#ifdef _WIN32
#include <Psapi.h>
#include <Windows.h>
#endif // _WIN32

#include "CmdLineParser.h"
#include "MarchingCubes.h"
#include "MemoryUsage.h"
#include "MultiGridOctreeData.h"
#include "Octree.h"
#include "PPolynomial.h"
#include "Ply.h"
#include "SparseMatrix.h"
#include "Time.h"

cmdLine<std::string> In("in");
cmdLine<std::string> Out("out");
cmdLine<std::string> VoxelGrid("voxel");
cmdLine<std::string> Xform("xForm");

#ifdef _WIN32
cmdLineReadable Performance("performance");
#endif
cmdLineReadable ShowResidual("showResidual");
cmdLineReadable NoComments("noComments");
cmdLineReadable PolygonMesh("polygonMesh");
cmdLineReadable Confidence("confidence");
cmdLineReadable NormalWeights("normalWeight");
cmdLineReadable NonManifold("nonManifold");
cmdLineReadable ASCII("ascii");
cmdLineReadable Density("density");
cmdLineReadable Verbose("verbose");

cmdLine<int> Depth("depth", 8);
cmdLine<int> SolverDivide("solverDivide", 8);
cmdLine<int> IsoDivide("isoDivide", 8);
cmdLine<int> KernelDepth("kernelDepth");
cmdLine<int> AdaptiveExponent("adaptiveExp", 1);
cmdLine<int> MinIters("minIters", 24);
cmdLine<int> FixedIters("iters", -1);
cmdLine<int> VoxelDepth("voxelDepth", -1);
#pragma message("[WARNING] Setting default min-depth to 5")
cmdLine<int> MinDepth("minDepth", 5);
cmdLine<int> MaxSolveDepth("maxSolveDepth" );
cmdLine<int> BoundaryType("boundary", 1);
cmdLine<int> Threads("threads", omp_get_num_procs());

cmdLine<float> SamplesPerNode("samplesPerNode", 1);
cmdLine<float> Scale("scale", 1.1);
cmdLine<float> SolverAccuracy("accuracy", 1e-3);
cmdLine<float> PointWeight("pointWeight", 4);

std::vector<cmdLineReadable*> params = {
	&In, &Depth, &Out, &Xform, &SolverDivide, &IsoDivide, &Scale, &Verbose, &SolverAccuracy, &NoComments,
	&KernelDepth, &SamplesPerNode, &Confidence, &NormalWeights, &NonManifold, &PolygonMesh, &ASCII,
	&ShowResidual, &MinIters, &FixedIters, &VoxelDepth, &PointWeight, &VoxelGrid, &Threads, &MinDepth,
	&MaxSolveDepth, &AdaptiveExponent, &BoundaryType, &Density,
#ifdef _WIN32
	&Performance
#endif
};

void ShowUsage(std::string const& executable) {
	printf( "Usage: %s\n" , executable.c_str() );
	printf( "\t --%s  <input points>\n" , In.name() );

	printf( "\t[--%s <ouput triangle mesh>]\n" , Out.name() );
	printf( "\t[--%s <ouput voxel grid>]\n" , VoxelGrid.name() );

	printf( "\t[--%s <maximum reconstruction depth>=%d]\n" , Depth.name() , Depth.value() );
	printf( "\t\t Running at depth d corresponds to solving on a 2^d x 2^d x 2^d\n" );
	printf( "\t\t voxel grid.\n" );

	printf( "\t[--%s <depth at which to extract the voxel grid>=<%s>]\n" , VoxelDepth.name() , Depth.name() );

	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name() , Scale.value() );
	printf( "\t\t Specifies the factor of the bounding cube that the input\n" );
	printf( "\t\t samples should fit into.\n" );

	printf( "\t[--%s <subdivision depth>=%d]\n" , SolverDivide.name() , SolverDivide.value() );
	printf( "\t\t The depth at which a block Gauss-Seidel solver is used\n");
	printf( "\t\t to solve the Laplacian.\n");

	printf( "\t[--%s <minimum number of samples per node>=%f]\n" , SamplesPerNode.name(), SamplesPerNode.value() );
	printf( "\t\t This parameter specifies the minimum number of points that\n" );
	printf( "\t\t should fall within an octree node.\n" );

	printf( "\t[--%s <num threads>=%d]\n" , Threads.name() , Threads.value() );
	printf( "\t\t This parameter specifies the number of threads across which\n" );
	printf( "\t\t the solver should be parallelizeds.\n" );

	printf( "\t[--%s]\n" , Confidence.name() );
	printf( "\t\t If this flag is enabled, the size of a sample's normals is\n" );
	printf( "\t\t used as a confidence value, affecting the sample's\n" );
	printf( "\t\t constribution to the reconstruction process.\n" );

	printf( "\t[--%s]\n" , NormalWeights.name() );
	printf( "\t\t If this flag is enabled, the size of a sample's normals is\n" );
	printf( "\t\t used as to modulate the interpolation weight.\n" );

	printf( "\t[--%s]\n" , NonManifold.name() );
	printf( "\t\t If this flag is enabled, the isosurface extraction does not add\n" );
	printf( "\t\t a planar polygon's barycenter in order to ensure that the output\n" );
	printf( "\t\t mesh is manifold.\n" );

	printf( "\t[--%s]\n" , PolygonMesh.name());
	printf( "\t\t If this flag is enabled, the isosurface extraction returns polygons\n" );
	printf( "\t\t rather than triangles.\n" );

	printf( "\t[--%s <interpolation weight>=%f]\n" , PointWeight.name() , PointWeight.value() );
	printf( "\t\t This value specifies the weight that point interpolation constraints are\n" );
	printf( "\t\t given when defining the (screened) Poisson system.\n" );

	printf( "\t[--%s <minimum depth>=%d]\n" , MinDepth.name() , MinDepth.value() );
	printf( "\t\t This flag specifies the minimum depth at which the octree is to be adaptive.\n" );

	printf( "\t[--%s <solver accuracy>=%g]\n" , SolverAccuracy.name() , SolverAccuracy.value ());
	printf( "\t[--%s <minimum number of solver iterations>=%d]\n" , MinIters.name() , MinIters.value() );

	printf( "\t[--%s <adaptive weighting exponent>=%d]\n", AdaptiveExponent.name() , AdaptiveExponent.value() );
	printf( "\t\t This flag specifies the exponent scale for the adaptive weighting.\n" );

#ifdef _WIN32
	printf( "\t[--%s]\n" , Performance.name() );
	printf( "\t\t If this flag is enabled, the running time and peak memory usage\n" );
	printf( "\t\t is output after the reconstruction.\n" );
#endif // _WIN32
	printf( "\t[--%s]\n" , Density.name() );
	printf( "\t[--%s]\n" , ASCII.name() );
	printf( "\t\t If this flag is enabled, the output file is written out in ASCII format.\n" );
	printf( "\t[--%s]\n" , NoComments.name() );
	printf( "\t\t If this flag is enabled, the output file will not include comments.\n" );
	printf( "\t[--%s]\n" , Verbose.name() );
	printf( "\t\t If this flag is enabled, the progress of the reconstructor will be output to STDOUT.\n" );
}

int ValidateFlags(std::string const& executable) {
	DumpOutput::instance().setEchoStdout(Verbose.set());
	DumpOutput::instance().setNoComments(NoComments.set());

	if(!In.set()) {
		ShowUsage(executable);
		return EXIT_FAILURE;
	}

	if(!MaxSolveDepth.set()) MaxSolveDepth.value() = Depth.value();

	if(SolverDivide.value() < MinDepth.value()) {
		std::cerr << "[WARNING] " << SolverDivide.name() << " must be at least as large as " <<
			MinDepth.name() << ": " << SolverDivide.value() << " >= " << MinDepth.value() << std::endl;
		SolverDivide.value() = MinDepth.value();
	}

	if(IsoDivide.value() < MinDepth.value()) {
		std::cerr << "[WARNING] " << IsoDivide.name() << " must be at least as large as " <<
			MinDepth.name() << ": " << IsoDivide.value() << " >= " << MinDepth.value() << std::endl;
		IsoDivide.value() = MinDepth.value();
	}

	if(!KernelDepth.set())
		KernelDepth.value() = Depth.value() - 2;

	if(KernelDepth.value() > Depth.value()) {
		std::cerr << "[ERROR] " << KernelDepth.name() << " can't be greater than " <<
			Depth.name() << ": " << KernelDepth.value() << " <= " << Depth.value();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

template<int Degree, class Real, class Vertex, bool OutputDensity>
int Execute() {
	DumpOutput::instance()("Running Screened Poisson Reconstruction (Version 5.71)\n");
	for(size_t i = 0; i != params.size(); ++i)
		if(params[i]->set())
			DumpOutput::instance()("\t--%s %s\n", params[i]->name(), params[i]->toString().c_str());

	XForm<Real, 4> xForm;
	if(Xform.set()) {
		std::ifstream file(Xform.value());
		if(!file) {
			std::cerr << "[WARNING] Could not read x-form from: " << Xform.value() << std::endl;
			xForm = XForm<Real, 4>::Identity();
		} else {
			for(int i = 0; i != 4; ++i)
				for(int j = 0; j != 4; ++j)
					file >> xForm(i, j);
		}
	}
	else xForm = XForm<Real, 4>::Identity();

	OctNode<TreeNodeData<OutputDensity>, Real>::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);

	double tt = Time();

	Octree<Degree, OutputDensity> tree(Threads.value(), Depth.value(), BoundaryType.value());

	double t = Time();
	tree.resetMaxMemoryUsage();
	int pointCount = tree.setTree(In.value().c_str(), Depth.value(), MinDepth.value(), KernelDepth.value(),
			SamplesPerNode.value(), Scale.value(), Confidence.set(), NormalWeights.set(), PointWeight.value(),
			AdaptiveExponent.value(), xForm);
	tree.ClipTree();
	tree.finalize(IsoDivide.value());

	DumpOutput::instance()("#             Tree set in: %9.1f (s), %9.1f (MB)\n", Time() - t,
			tree.maxMemoryUsage());
	DumpOutput::instance()("#               Input Points: %d\n", pointCount);
	DumpOutput::instance()("#               Leaves/Nodes: %lld/%lld\n", tree.tree().leaves(),
			tree.tree().nodes());
	DumpOutput::instance()("#               Memory Usage: %.3f MB\n",
			float(MemoryInfo::Usage()) / (1 << 20));

	double maxMemoryUsage = tree.maxMemoryUsage();
	t = Time();
	tree.resetMaxMemoryUsage();
	tree.SetLaplacianConstraints();
	DumpOutput::instance()("#      Constraints set in: %9.1f (s), %9.1f (MB)\n", Time() - t,
			tree.maxMemoryUsage());
	DumpOutput::instance()("#               Memory Usage: %.3f MB\n",
			float(MemoryInfo::Usage()) / (1 << 20));
	maxMemoryUsage = std::max(maxMemoryUsage, tree.maxMemoryUsage());

	t = Time();
	tree.resetMaxMemoryUsage();
	tree.LaplacianMatrixIteration(SolverDivide.value(), ShowResidual.set(), MinIters.value(),
			SolverAccuracy.value(), MaxSolveDepth.value(), FixedIters.value());
	DumpOutput::instance()("# Linear system solved in: %9.1f (s), %9.1f (MB)\n", Time() - t,
			tree.maxMemoryUsage());
	DumpOutput::instance()("#            Memory Usage: %.3f MB\n", float(MemoryInfo::Usage()) / (1 << 20));
	maxMemoryUsage = std::max(maxMemoryUsage, tree.maxMemoryUsage());

	t = Time();
	Real isoValue = tree.GetIsoValue();
	DumpOutput::instance()("#          Got average in: %f\n", Time() - t);
	DumpOutput::instance()("#               Iso-Value: %e\n", isoValue);

	if(VoxelGrid.set()) {
		double t = Time();
		std::ofstream file(VoxelGrid.value(), std::ofstream::out | std::ofstream::binary);
		if(!file) std::cerr << "Failed to open voxel file for writing: " << VoxelGrid.value() << std::endl;
		else {
			int res;
			Pointer(Real) values = tree.GetSolutionGrid(res, isoValue, VoxelDepth.value());
			file.write(reinterpret_cast<char*>(&res), sizeof(res));
			for(int i = 0; i != res * res * res; ++i) {
				float v = (float)values[i];
				file.write(reinterpret_cast<char*>(&v), sizeof(v));
			}
			DeletePointer(values);
		}
		DumpOutput::instance()("#       Got voxel grid in: %f\n" , Time()-t );
	}

	if(Out.set()) {
		t = Time();
		CoredFileMeshData<Vertex> mesh;
		tree.resetMaxMemoryUsage();
		tree.GetMCIsoTriangles(isoValue, IsoDivide.value(), &mesh, 0, 1, !NonManifold.set(),
				PolygonMesh.set());
		if(PolygonMesh.set())
			DumpOutput::instance()("#         Got polygons in: %9.1f (s), %9.1f (MB)\n", Time() - t,
					tree.maxMemoryUsage());
		else
			DumpOutput::instance()("#        Got triangles in: %9.1f (s), %9.1f (MB)\n", Time() - t,
					tree.maxMemoryUsage());
		maxMemoryUsage = std::max(maxMemoryUsage, tree.maxMemoryUsage());
		DumpOutput::instance()("#             Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - tt,
				maxMemoryUsage);

		PlyWritePolygons(Out.value().c_str(), &mesh, ASCII.set() ? PLY_ASCII : PLY_BINARY_NATIVE,
				DumpOutput::instance().strings(), xForm.inverse());
	}

	return EXIT_SUCCESS;
}

#ifdef _WIN32
inline double to_seconds( const FILETIME& ft )
{
	const double low_to_sec=100e-9; // 100 nanoseconds
	const double high_to_sec=low_to_sec*4294967296.0;
	return ft.dwLowDateTime*low_to_sec+ft.dwHighDateTime*high_to_sec;
}
#endif // _WIN32

int main(int argc, char** argv) {
	cmdLineParse(argc - 1, argv + 1, params);
	int ret;
	if((ret = ValidateFlags(argv[0]))) return ret;
	ret = Density.set() ? Execute<2, Real, PlyValueVertex<Real>, true>() :
		Execute<2, Real, PlyVertex<Real>, false>();
#ifdef _WIN32
	if( Performance.set() )
	{
		HANDLE cur_thread=GetCurrentThread();
		FILETIME tcreat, texit, tkernel, tuser;
		if( GetThreadTimes( cur_thread , &tcreat , &texit , &tkernel , &tuser ) )
			printf( "Time (Wall/User/Kernel): %.2f / %.2f / %.2f\n" , Time()-t , to_seconds( tuser ) , to_seconds( tkernel ) );
		else printf( "Time: %.2f\n" , Time()-t );
		HANDLE h = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		if( GetProcessMemoryInfo( h , &pmc , sizeof(pmc) ) ) printf( "Peak Memory (MB): %d\n" , pmc.PeakWorkingSetSize>>20 );
	}
#endif // _WIN32
	return ret;
}
