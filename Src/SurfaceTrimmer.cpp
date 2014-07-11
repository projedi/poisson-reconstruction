/*
Copyright (c) 2013, Michael Kazhdan
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

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <iostream>

#ifndef NO_OMP
#include <omp.h>
#endif

#include "CmdLineParser.h"
#include "DumpOutput.h"
#include "HashMap.h"
#include "Geometry.h"
#include "Ply.h"
#include "MAT.h"
#include "Time.h"

#define FOR_RELEASE 1

cmdLine<std::string> In("in");
cmdLine<std::string> Out("out");
cmdLine<int> Smooth("smooth", 5);
cmdLine<float> Trim("trim");
cmdLine<float> IslandAreaRatio("aRatio", 0.001f);
cmdLine<std::pair<float, float> > ColorRange("color");
cmdLineReadable PolygonMesh("polygonMesh");

std::vector<cmdLineReadable*> params;

void BuildParams() {
	cmdLineReadable* params_array[] = {
		&In , &Out , &Trim , &PolygonMesh , &ColorRange , &Smooth , &IslandAreaRatio, nullptr
	};
	for(cmdLineReadable** p = params_array; *p; ++p)
		params.push_back(*p);
}

void ShowUsage(std::string const& executable) {
	printf( "Usage: %s\n" , executable.c_str() );
	printf( "\t --%s <input polygon mesh>\n" , In.name() );
	printf( "\t[--%s <ouput polygon mesh>]\n" , Out.name() );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , Smooth.name() , Smooth.value() );
	printf( "\t[--%s <trimming value>]\n" , Trim.name() );
	printf( "\t[--%s <relative area of islands>=%f]\n" , IslandAreaRatio.name() , IslandAreaRatio.value() );
	printf( "\t[--%s]\n" , PolygonMesh.name() );
#if !FOR_RELEASE
	printf( "\t[--%s <color range>]\n" , ColorRange.name() );
#endif // !FOR_RELEASE
}

long long EdgeKey(long long key1, long long key2) {
	return key1 <= key2 ? (key1 << 32) | key2 : EdgeKey(key2, key1);
}

template<class Real>
PlyValueVertex<Real> InterpolateVertices(PlyValueVertex<Real> const& v1, PlyValueVertex<Real> const& v2,
		float value) {
	if(v1.value == v2.value) return (v1 + v2) / (Real)2;

	Real dx = (v1.value - value) / (v1.value - v2.value);
	PlyValueVertex<Real> v;
	for(int i = 0; i != 3; ++i) v.point.coords[i] = v1.point.coords[i] * (1 - dx) + v2.point.coords[i] * dx;
	v.value = v1.value * (1 - dx) + v2.value * dx;
	return v;
}

template<class Real>
std::vector<PlyColorVertex<Real> > ColorVertices(std::vector<PlyValueVertex<Real> > const& inVertices,
		float min, float max) {
	std::vector<PlyColorVertex<Real> > outVertices(inVertices.size());
	for(size_t i = 0; i != inVertices.size(); ++i) {
		outVertices[i].point = inVertices[i].point;
		float temp = (inVertices[i].value - min) / (max - min);
		temp = std::max(0.f, std::min(1.f, temp));
		temp *= 255;
		outVertices[i].color[0] = outVertices[i].color[1] = outVertices[i].color[2] = (int)std::lround(temp);
	}
	return outVertices;
}

template<class Real>
void SmoothValues(std::vector<PlyValueVertex<Real> >& vertices,
		std::vector<std::vector<int> > const& polygons) {
	std::vector<int> count(vertices.size());
	std::vector<Real> sums(vertices.size(), 0);
	for(size_t i = 0; i != polygons.size(); ++i) {
		size_t sz = polygons[i].size();
		for(size_t j = 0; j != sz; ++j) {
			int v1 = polygons[i][j];
			int v2 = polygons[i][(j + 1) % sz];
			++count[v1];
			++count[v2];
			sums[v1] += vertices[v2].value;
			sums[v2] += vertices[v1].value;
		}
	}
	for(size_t i = 0; i != vertices.size(); ++i)
		vertices[i].value = (sums[i] + vertices[i].value) / (count[i] + 1);
}

template<class Real>
void SplitPolygon(std::vector<int> const& polygon, std::vector<PlyValueVertex<Real> >& vertices, 
		std::vector<std::vector<int> >& ltPolygons, std::vector<std::vector<int> >& gtPolygons,
		std::vector<bool>& ltFlags, std::vector<bool>& gtFlags,
		HashMap<long long, int>& vertexTable, Real trimValue) {
	int sz = polygon.size();
	std::vector<bool> gt(sz);
	int gtCount = 0;
	for(int j = 0; j != sz; ++j) {
		gt[j] = vertices[polygon[j]].value > trimValue ;
		if(gt[j]) ++gtCount;
	}
	if(gtCount == sz) {
		gtPolygons.push_back(polygon);
		gtFlags.push_back(false);
	} else if(gtCount == 0) {
		ltPolygons.push_back(polygon);
		ltFlags.push_back(false);
	} else {
		int start;
		for(start = 0; start != sz; ++start) if(gt[start] && !gt[(start + sz - 1) % sz]) break;

		bool gtFlag = true;
		std::vector<int> poly;

		// Add the initial vertex
		{
			int v1 = polygon[(start + sz - 1) % sz];
			int v2 = polygon[start];
			int vIdx;
			HashMap<long long, int>::iterator iter = vertexTable.find(EdgeKey(v1, v2));
			if(iter == vertexTable.end()) {
				vertexTable[EdgeKey(v1, v2)] = vIdx = vertices.size();
				vertices.push_back(InterpolateVertices(vertices[v1], vertices[v2], trimValue));
			} else vIdx = iter->second;
			poly.push_back(vIdx);
		}

		for(int _j = 0; _j <= sz; ++_j) {
			int j1 = (_j + start + sz - 1) % sz;
			int j2 = (_j + start) % sz;
			int v1 = polygon[j1];
			int v2 = polygon[j2];
			if(gt[j2] == gtFlag) poly.push_back(v2);
			else {
				int vIdx;
				HashMap<long long, int>::iterator iter = vertexTable.find(EdgeKey(v1, v2));
				if(iter == vertexTable.end()) {
					vertexTable[EdgeKey(v1, v2)] = vIdx = vertices.size();
					vertices.push_back(InterpolateVertices(vertices[v1], vertices[v2], trimValue));
				} else vIdx = iter->second;
				poly.push_back(vIdx);
				if(gtFlag) {
					gtPolygons.push_back(poly);
					ltFlags.push_back(true);
				} else {
					ltPolygons.push_back(poly);
					gtFlags.push_back(true);
				}
				poly.clear();
				poly.push_back(vIdx);
				poly.push_back(v2);
				gtFlag = !gtFlag;
			}
		}
	}
}

template<class Real>
std::vector<std::vector<int> > Triangulate(std::vector<PlyValueVertex<Real> > const& vertices,
		std::vector<std::vector<int> > const& polygons) {
	std::vector<std::vector<int> > triangles;
	for(size_t i = 0; i != polygons.size(); ++i) {
		if(polygons.size() > 3) {
			MinimalAreaTriangulation<Real> mat;
			std::vector<Point3D<Real> > _vertices(polygons[i].size());
			std::vector<TriangleIndex> _triangles;
			for(size_t j = 0; j != polygons[i].size(); ++j)
				_vertices[j] = vertices[polygons[i][j]].point;
			mat.GetTriangulation(_vertices, _triangles);

			// Add the triangles to the mesh
			size_t idx = triangles.size();
			triangles.resize(idx + _triangles.size());
			for(size_t j = 0; j != _triangles.size(); ++j) {
				triangles[idx + j].resize(3);
				for(int k = 0; k != 3; ++k) triangles[idx + j][k] = polygons[i][_triangles[j].idx[k]];
			}
		} else if(polygons[i].size() == 3) triangles.push_back(polygons[i]);
	}
	return triangles;
}

template<class Vertex>
void RemoveHangingVertices(std::vector<Vertex>& vertices, std::vector<std::vector<int> >& polygons) {
	HashMap<int, int> vMap;
	std::vector<bool> vertexFlags(vertices.size(), false);
	for(size_t i = 0; i != polygons.size(); ++i)
		for(size_t j = 0; j != polygons[i].size(); ++j)
			vertexFlags[polygons[i][j]] = true;
	int vCount = 0;
	for(size_t i = 0; i != vertices.size(); ++i)
		if(vertexFlags[i]) vMap[i] = vCount++;
	for(size_t i = 0; i != polygons.size(); ++i)
		for(size_t j = 0; j != polygons[i].size(); ++j)
			polygons[i][j] = vMap[polygons[i][j]];

	std::vector<Vertex> _vertices(vCount);
	for(size_t i = 0; i != vertices.size(); ++i)
		if(vertexFlags[i]) _vertices[vMap[i]] = vertices[i];
	vertices = _vertices;
}

std::vector<std::vector<int> > SetConnectedComponents(std::vector<std::vector<int> > const& polygons) {
	std::vector<std::vector<int> > components;
	std::vector<int> polygonRoots(polygons.size());
	for(size_t i = 0; i != polygons.size(); ++i)
		polygonRoots[i] = i;
	HashMap<long long, int> edgeTable;
	for(size_t i = 0; i != polygons.size(); ++i) {
		int sz = polygons[i].size();
		for(int j = 0; j != sz; ++j) {
			long long eKey = EdgeKey(polygons[i][j], polygons[i][(j + 1) % sz]);
			HashMap<long long, int>::iterator iter = edgeTable.find(eKey);
			if(iter == edgeTable.end()) edgeTable[eKey] = i;
			else {
				int p = iter->second;
				while(polygonRoots[p] != p) {
					int temp = polygonRoots[p];
					polygonRoots[p] = i;
					p = temp;
				}
				polygonRoots[p] = i;
			}
		}
	}
	for(size_t i = 0; i != polygonRoots.size(); ++i) {
		int p = i;
		while(polygonRoots[p] != p) p = polygonRoots[p];
		int root = p;
		p = i;
		while(polygonRoots[p] != p) {
			int temp = polygonRoots[p];
			polygonRoots[p] = root;
			p = temp;
		}
	}
	int cCount = 0;
	HashMap<int, int> vMap;
	for(size_t i = 0; i != polygonRoots.size(); ++i) if(polygonRoots[i] == (int)i) vMap[i] = cCount++;
	components.resize(cCount);
	for(size_t i = 0; i != polygonRoots.size(); ++i) components[vMap[polygonRoots[i]]].push_back(i);
	return components;
}

template<class Real>
double PolygonArea(std::vector<PlyValueVertex<Real> > const& vertices, std::vector<int> const& polygon) {
	if(polygon.size() < 3) return 0;
	else if(polygon.size() == 3)
		return TriangleArea(vertices[polygon[0]].point, vertices[polygon[1]].point,
				vertices[polygon[2]].point);
	else {
		Point3D<Real> center;
		for(size_t i = 0; i != polygon.size(); ++i)
			center += vertices[polygon[i]].point;
		center /= (Real)polygon.size();
		double area = 0;
		for(size_t i = 0; i != polygon.size(); ++i)
			area += TriangleArea(center, vertices[polygon[i]].point, vertices[polygon[(i + 1) % polygon.size()]].point);
		return area;
	}
}

int main(int argc, char** argv) {
	BuildParams();
	cmdLineParse(argc - 1, argv + 1, params , false);

#if FOR_RELEASE
	if(!In.set() || !Trim.set()) {
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
#else // !FOR_RELEASE
	if(!In.set()) {
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
#endif // FOR_RELEASE
	std::vector<PlyValueVertex<float> > vertices;
	std::vector<std::vector<int> > polygons;
	int ft;
	std::vector<std::string> comments;
	bool readFlags[PlyValueVertex<float>::Components];
	PlyReadPolygons(In.value(), vertices, polygons, ft, comments, readFlags);
	DumpOutput::instance().resetStrings(comments);
	if(!readFlags[3]) {
		std::cerr << "[ERROR] vertices do not have value flag" << std::endl;
		return EXIT_FAILURE;
	}
	for(int i = 0; i != Smooth.value(); ++i) SmoothValues(vertices, polygons);
	float min = vertices[0].value;
	float max = vertices[0].value;
	for(size_t i = 0; i != vertices.size(); ++i) {
		min = std::min(min, vertices[i].value);
		max = std::max(max, vertices[i].value);
	}
	std::cout << "Value Range: [" << min << ", " << max << "]" << std::endl;

	if(Trim.set()) {
		DumpOutput::instance()("Running Surface Trimmer (V5)");
		for(auto p: params)
			if(p->set())
				DumpOutput::instance()("\t--%s %s\n", p->name(), p->toString().c_str());

		HashMap<long long, int> vertexTable;
		std::vector<std::vector<int> > ltPolygons;
		std::vector<std::vector<int> > gtPolygons;
		std::vector<bool> ltFlags;
		std::vector<bool> gtFlags;

		double t = Time();
		for(size_t i = 0; i != polygons.size(); ++i)
			SplitPolygon(polygons[i], vertices, ltPolygons, gtPolygons, ltFlags, gtFlags,
					vertexTable, Trim.value());
		if(IslandAreaRatio.value() > 0) {
			std::vector<std::vector<int> > _gtPolygons;
			std::vector<std::vector<int> > ltComponents = SetConnectedComponents(ltPolygons);
			std::vector<std::vector<int> > gtComponents = SetConnectedComponents(gtPolygons);
			std::vector<double> ltAreas(ltComponents.size(), 0);
			std::vector<double> gtAreas(gtComponents.size(), 0);
			std::vector<bool> ltComponentFlags(ltComponents.size(), false);
			std::vector<bool> gtComponentFlags(gtComponents.size(), false);
			double area = 0;
			for(size_t i = 0; i != ltComponents.size(); ++i) {
				for(size_t j = 0; j != ltComponents[i].size(); ++j) {
					ltAreas[i] += PolygonArea(vertices, ltPolygons[ltComponents[i][j]]);
					ltComponentFlags[i] = ltComponentFlags[i] || ltFlags[ltComponents[i][j]];
				}
				area += ltAreas[i];
			}
			for(size_t i = 0; i != gtComponents.size(); ++i) {
				for(size_t j = 0; j != gtComponents[i].size(); ++j) {
					gtAreas[i] += PolygonArea(vertices, gtPolygons[gtComponents[i][j]]);
					gtComponentFlags[i] = gtComponentFlags[i] || gtFlags[gtComponents[i][j]];
				}
				area += gtAreas[i];
			}
			for(size_t i = 0; i != ltComponents.size(); ++i) {
				if(ltAreas[i] < area * IslandAreaRatio.value() && ltComponentFlags[i]) {
					for(size_t j = 0; j != ltComponents[i].size(); ++j)
						_gtPolygons.push_back(ltPolygons[ltComponents[i][j]]);
				}
			}
			for(size_t i = 0; i != gtComponents.size(); ++i) {
				if(gtAreas[i] >= area * IslandAreaRatio.value() && gtComponentFlags[i]) {
					for(size_t j = 0; j != gtComponents[i].size(); ++j)
						_gtPolygons.push_back(gtPolygons[gtComponents[i][j]]);
				}
			}
			gtPolygons = _gtPolygons;
		}
		std::vector<std::vector<int> > polys =
			PolygonMesh.set() ? gtPolygons : Triangulate(vertices, gtPolygons);

		RemoveHangingVertices(vertices, polys);
		DumpOutput::instance()("#Trimmed In: %9.1f (s)", Time() - t);
		if(Out.set())
			PlyWritePolygons(Out.value(), vertices, polys, ft, DumpOutput::instance().strings());
	} else {
		if(ColorRange.set()) {
			min = ColorRange.value().first;
			max = ColorRange.value().second;
		}
		std::vector<PlyColorVertex<float> > outVertices = ColorVertices(vertices, min, max);
		if(Out.set())
			PlyWritePolygons(Out.value(), outVertices, polygons, ft, DumpOutput::instance().strings());
	}
	return EXIT_SUCCESS;
}
