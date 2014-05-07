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

#include <cstdio>
#include <cstring>
#include <iostream>

#include "Geometry.h"

// TODO: fread, fwrite and the kernel already have an underlying buffering.
// TODO: So remove that with passion.

BufferedReadWriteFile::BufferedReadWriteFile():
	buffer_index_(0),
	buffer_size_(1 << 20) {
#ifdef _WIN32
	tmpfile_s(&fp_);
#else
	fp_ = tmpfile();
#endif
	if(errno) {
		perror("[ERROR] Failed to create temporary file\n");
		exit(1);
	}
	buffer_ = (char*)malloc(buffer_size_);
}

BufferedReadWriteFile::~BufferedReadWriteFile() {
	free(buffer_);
	fclose(fp_);
}

void BufferedReadWriteFile::reset() {
	if(buffer_index_)
		fwrite(buffer_, 1, buffer_index_, fp_);
	buffer_index_ = 0;
	fseek(fp_, 0, SEEK_SET);
	buffer_index_ = 0;
	buffer_size_ = fread(buffer_, 1, buffer_size_, fp_);
}

bool BufferedReadWriteFile::write(void const* data, size_t size) {
	if(!size) return true;
	char* _data = (char*)data;
	size_t sz = buffer_size_ - buffer_index_;
	while(sz <= size) {
		memcpy(buffer_ + buffer_index_, _data, sz);
		fwrite(buffer_, 1, buffer_size_, fp_);
		_data += sz;
		size -= sz;
		buffer_index_ = 0;
		sz = buffer_size_;
	}
	if(size) {
		memcpy(buffer_ + buffer_index_, _data, size);
		buffer_index_ += size;
	}
	return true;
}

bool BufferedReadWriteFile::read(void* data, size_t size) {
	if(!size) return true;
	char *_data = (char*)data;
	size_t sz = buffer_size_ - buffer_index_;
	while(sz <= size) {
		if(size && !buffer_size_) return false;
		memcpy(_data, buffer_ + buffer_index_, sz);
		buffer_size_ = fread(buffer_, 1, buffer_size_, fp_);
		_data += sz;
		size -= sz;
		buffer_index_ = 0;
		if(!size) return true;
		sz = buffer_size_;
	}
	if(size) {
		if(!buffer_size_) return false;
		memcpy(_data, buffer_ + buffer_index_, size);
		buffer_index_ += size;
	}
	return true;
}
