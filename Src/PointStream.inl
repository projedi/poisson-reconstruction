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

template<class Real>
ASCIIPointStream<Real>::ASCIIPointStream(std::string const& filename): file_(filename) {
	if(!file_) {
		std::cerr << "Failed to open file for reading: " << filename << std::endl;
		std::exit(1);
	}
}

template<class Real>
bool ASCIIPointStream<Real>::nextPoint(Point3D<Real>& p, Point3D<Real>& n) {
	return (file_ >> p[0] >> p[1] >> p[2] >> n[0] >> n[1] >> n[2]).good();
}

template<class Real>
BinaryPointStream<Real>::BinaryPointStream(std::string const& filename):
	file_(filename, std::ios_base::binary) {
	if(!file_) {
		std::cerr << "Failed to open file for reading: " << filename << std::endl;
		std::exit(1);
	}
}


template<class Real>
bool BinaryPointStream<Real>::nextPoint(Point3D<Real>& p, Point3D<Real>& n) {
	Real buf[6];
	file_.read((char*)buf, 6 * sizeof(Real));
	p[0] = buf[0];
	p[1] = buf[1];
	p[2] = buf[2];
	n[0] = buf[3];
	n[1] = buf[4];
	n[2] = buf[5];
	return file_.good();
}

template<class Real>
void PLYPointStream<Real>::reset() {
	int fileType;
	float version;
	PlyProperty** plist;
	if(_ply) _free();
	_ply = ply_open_for_reading(_fileName.c_str(), &_nr_elems, &_elist, &fileType, &version);
	if(!_ply) {
		std::cerr << "[ERROR] Failed to open ply file for reading: " << _fileName << std::endl;
		std::exit(1);
	}
	bool foundVertices = false;
	for(int i = 0; i != _nr_elems; ++i) {
		int num_elems;
		int nr_props;
		char* elem_name = _elist[i];
		plist = ply_get_element_description(_ply , elem_name , &num_elems , &nr_props);
		if(!plist) {
			std::cerr << "[ERROR] Failed to get element description: " << elem_name << std::endl;
			std::exit(1);
		}

		if(equal_strings("vertex" , elem_name)) {
			foundVertices = true;
			_pCount = num_elems;
			_pIdx = 0;
			for(int i = 0; i != PlyOrientedVertex<Real>::Components; ++i) 
				if(!ply_get_property(_ply, elem_name, &(PlyOrientedVertex<Real>::Properties[i]))) {
					std::cerr << "[ERROR] Failed to find property in ply file: " <<
						PlyOrientedVertex< Real >::Properties[i].name << std::endl;
					std::exit(1);
				}
		}
		for(int j = 0; j != nr_props; ++j) {
			free(const_cast<char*>(plist[j]->name));
			free(plist[j]);
		}
		free(plist);
		if(foundVertices) break;
	}
	if(!foundVertices) {
		std::cerr << "[ERROR] Could not find vertices in ply file" << std::endl;
		std::exit(1);
	}
}

template<class Real>
void PLYPointStream<Real>::_free() {
	if(_ply) {
		ply_close(_ply);
		_ply = NULL;
	}
	if(_elist) {
		for(int i = 0; i != _nr_elems; ++i) free(_elist[i]);
		free(_elist);
	}
}

template<class Real>
bool PLYPointStream<Real>::nextPoint(Point3D<Real>& p, Point3D<Real>& n) {
	if(_pIdx >= _pCount) return false;
	PlyOrientedVertex<Real> op;
	ply_get_element(_ply, (void *)&op);
	p = op.point;
	n = op.normal;
	++_pIdx;
	return true;
}

bool strcaseequal(std::string const& s1, std::string const& s2) {
#ifdef WIN32
	int res = _stricmp(s1.c_str(), s2.c_str());
#else
	int res = strcasecmp(s1.c_str(), s2.c_str());
#endif
	return !res;
}

template<class Real>
PointStream<Real>* PointStream<Real>::open(std::string const& filename) {
	size_t last_dot = filename.find_last_of('.');
	std::string ext = last_dot == std::string::npos ? "" : filename.substr(last_dot + 1);
	if(strcaseequal(ext, "bnpts")) return new BinaryPointStream<Real>(filename);
	if(strcaseequal(ext, "ply")) return new PLYPointStream<Real>(filename);
	else return new ASCIIPointStream<Real>(filename);
}
