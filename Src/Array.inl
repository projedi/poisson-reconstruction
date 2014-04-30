/*
Copyright (c) 2011, Michael Kazhdan and Ming Chuang
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

#define FULL_ARRAY_DEBUG 0 // Note that this is not thread-safe

#include <cstdio>
#include <cstddef>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef _WIN64
#define ASSERT(x) { if(!(x)) __debugbreak(); }
#elif defined(_WIN32)
#define ASSERT(x) { if(!(x)) _asm{ int 0x03 } }
#else
#define ASSERT(x) { if(!(x)) exit(0); }
#endif

#if FULL_ARRAY_DEBUG
struct DebugMemoryInfo {
	void const* address;
	char name[512];
};
static std::vector<DebugMemoryInfo> memoryInfo;
#endif

template<class C>
class Array {
	void _assertBounds(long long idx) const {
		if(idx < min || idx >= max) {
			std::cerr << "Array index out-of-bounds: " << min << " <= " << idx
				<< " < " << max << std::endl;
			ASSERT(0);
			exit(0);
		}
	}
protected:
	C *data , *_data;
	long long min;
	long long max;
#if FULL_ARRAY_DEBUG
	static void _AddMemoryInfo(void const* ptr, char const* name) {
		size_t sz = memoryInfo.size();
		memoryInfo.resize(sz + 1);
		memoryInfo[sz].address = ptr;
		if(name) strcpy(memoryInfo[sz].name, name);
		else memoryInfo[sz].name[0] = 0;
	}

	static void _RemoveMemoryInfo(void const* ptr) {
		size_t idx;
		for(idx = 0; idx != memoryInfo.size(); ++idx)
			if(memoryInfo[idx].address == ptr) break;
		if(idx == memoryInfo.size()) {
			std::cerr << "Could not find memory address table" << std::endl;
			ASSERT(0);
		} else {
			memoryInfo[idx] = memoryInfo[memoryInfo.size() - 1];
			memoryInfo.pop_back();
		}
	}
#else
	static void _AddMemoryInfo(void const*, char const*) { }
	static void _RemoveMemoryInfo(void const*) { }
#endif // FULL_ARRAY_DEBUG

public:
	long long minimum() const { return min; }
	long long maximum() const { return max; }

	static Array New(size_t size , char const* name = nullptr) {
		Array a;
		a._data = a.data = new C[size];
		a.min = 0;
#pragma message( "[WARNING] Casting unsigned to signed" )
		a.max = (long long) size;
		_AddMemoryInfo(a._data, name );
		return a;
	}

	static Array Alloc(size_t size, bool clear, char const* name = nullptr) {
		Array a;
		a._data = a.data = (C*)malloc(size * sizeof(C));
		if(clear) memset(a.data, 0, size * sizeof(C));
		a.min = 0;
#pragma message( "[WARNING] Casting unsigned to signed" )
		a.max = (long long) size;
		_AddMemoryInfo(a._data, name);
		return a;
	}

	Array(): data(nullptr), _data(nullptr), min(0), max(0) { }

	template<class D>
	Array(Array<D>& a): _data(nullptr) {
		if(!a) {
			data = nullptr;
			min = max = 0;
		} else {
			// [WARNING] Chaning szC and szD to size_t causes some really strange behavior.
			long long szC = sizeof(C);
			long long szD = sizeof(D);
			data = (C*)&a[0];
			min = (a.minimum() * szD) / szC;
			max = (a.maximum() * szD) / szC;
			if(min * szC != a.minimum() * szD || max * szC != a.maximum() * szD) {
				std::cerr << "Could not convert array [ " << a.minimum() << " , "
					<< a.maximum() << " ] * " << szD << " => [ " << min << " , "
					<< max << " * szC" << std::endl;
				ASSERT(0);
				exit(0);
			}
		}
	}

	bool operator==(Array<C> const& a) const { return data == a.data; }
	bool operator!=(Array<C> const& a) const { return !(*this == a); }
	bool operator==(C const* c) const { return data == c; }
	bool operator!=(C const* c) const { return !(*this == c); }

	C* operator->() {
		return const_cast<C*>(const_cast<Array const*>(this)->operator->());
	}

	C const* operator->() const {
		_assertBounds(0);
		return data;
	}

	C& operator[](long long idx) {
		return const_cast<C&>(const_cast<Array const*>(this)->operator[](idx));
	}

	C const& operator[](long long idx) const {
		_assertBounds(idx);
		return data[idx];
	}

	Array operator+(long long idx) const {
		Array a;
		a._data = _data;
		a.data = data + idx;
		a.min = min - idx;
		a.max = max - idx;
		return a;
	}

	Array& operator+=(long long idx) {
		min -= idx;
		max -= idx;
		data += idx;
		return *this;
	}

	Array& operator++() { return *this += 1; }

	Array operator-(long long idx) const { return *this + (-idx); }

	Array& operator-=(long long idx) { return *this += -idx; }

	Array& operator--() { return *this -= 1; }

	long long operator-(Array const& a) const { return data - a.data; }

	void Free() {
		if(_data) {
			free(_data);
			_RemoveMemoryInfo(_data);
		}
		*this = Array();
	}

	void Delete() {
		if(_data) {
			delete[] _data;
			_RemoveMemoryInfo(_data);
		}
		*this = Array();
	}

	C* pointer() { return const_cast<C*>(const_cast<Array const*>(this)->pointer()); }
	C const* pointer() const { return data; }

	bool operator!() const { return data == nullptr; }
	operator bool() const { return data != nullptr; }
};

template<class C>
class ConstArray {
	void _assertBounds(long long idx) const {
		if(idx < min || idx >= max) {
			std::cerr << "Array index out-of-bounds: " << min << " <= " << idx
				<< " < " << max << std::endl;
			ASSERT(0);
			exit(0);
		}
	}
protected:
	C const* data;
	long long min;
	long long max;
public:
	long long minimum() const { return min; }
	long long maximum() const { return max; }

	ConstArray(): data(nullptr), min(0), max(0) { }

	ConstArray(Array<C> const& a):
		data(a.pointer()), min(a.minimum()), max(a.maximum()) { }

	template<class D>
	ConstArray(Array<D> const& a) {
		// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
		long long szC = sizeof(C);
		long long szD = sizeof(D);
		data = (const C*)a.pointer();
		min = (a.minimum() * szD) / szC;
		max = (a.maximum() * szD) / szC;
		if(min * szC != a.minimum() * szD || max * szC != a.maximum() * szD) {
			std::cerr << "Could not convert array [ " << a.minimum() << " , "
				<< a.maximum() << " ] * " << szD << " => [ " << min << " , "
				<< max << " * szC" << std::endl
				<< a.minimum() << " " << a.minimum() * szD << " " << (a.minimum() * szD) / szC
				<< std::endl;
			ASSERT(0);
			exit(0);
		}
	}

	template<class D>
	ConstArray(ConstArray<D> const& a) {
		// [WARNING] Chaning szC and szD to size_t causes some really strange behavior.
		long long szC = sizeof(C);
		long long szD = sizeof(D);
		data = (const C*)a.pointer();
		min = (a.minimum() * szD) / szC;
		max = (a.maximum() * szD) / szC;
		if(min * szC != a.minimum() * szD || max * szC != a.maximum() * szD) {
			std::cerr << "Could not convert array [ " << a.minimum() << " , "
				<< a.maximum() << " ] * " << szD << " => [ " << min << " , "
				<< max << " * szC" << std::endl;
			ASSERT(0);
			exit(0);
		}
	}

	bool operator==(ConstArray<C> const& a) const { return data == a.data; }
	bool operator!=(ConstArray<C> const& a) const { return !(*this == a); }
	bool operator==(C const* c) const { return data == c; }
	bool operator!=(C const* c) const { return !(*this == c); }

	C const* operator->() {
		_assertBounds(0);
		return data;
	}

	C const& operator[](long long idx) const {
		_assertBounds(idx);
		return data[idx];
	}

	ConstArray operator+(long long idx) const {
		ConstArray a;
		a.data = data + idx;
		a.min = min - idx;
		a.max = max - idx;
		return a;
	}

	ConstArray& operator+=(long long idx) {
		min -= idx;
		max -= idx;
		data += idx;
		return *this;
	}

	ConstArray& operator++() { return *this += 1; }

	ConstArray operator-(long long idx) const { return *this + (-idx); }

	ConstArray& operator-=(long long idx) { return *this += (-idx); }

	ConstArray& operator--() { return *this -= 1; }

	long long operator-(ConstArray const& a) const { return data - a.data; }
	long long operator-(Array<C> const& a) const { return data - a.pointer(); }

	C const* pointer() const { return data; }
	bool operator !() { return data == nullptr; }
	operator bool() { return data != nullptr; }
};

#if FULL_ARRAY_DEBUG
inline void PrintMemoryInfo() {
	for(size_t i = 0; i != memoryInfo.size(); ++i)
		std::cout << i << "] " << memoryInfo[i].name << std::endl;
}
#endif // FULL_ARRAY_DEBUG

template<class C>
Array<C> memcpy(Array<C> destination, void const* source, size_t size) {
	if(size > destination.maximum() * sizeof(C)) {
		std::cerr << "Size of copy exceeds destination maximum: " << size
			<< " > " << destination.maximum() * sizeof(C) << std::endl;
		ASSERT(0);
		exit(0);
	}
	memcpy(&destination[0], source, size);
	return destination;
}

template<class C, class D>
Array<C> memcpy(Array<C> destination, Array<D> source, size_t size) {
	if(size > source.maximum() * sizeof(D)) {
		std::cerr << "Size of copy exceeds source maximum: " << size
			<< " > " << source.maximum() * sizeof(D) << std::endl;
		ASSERT(0);
		exit(0);
	}
	memcpy(destination, &source[0], size);
	return destination;
}

template<class C, class D>
Array<C> memcpy(Array<C> destination, ConstArray<D> source, size_t size) {
	if(size > source.maximum() * sizeof(D)) {
		std::cerr << "Size of copy exceeds source maximum: " << size
			<< " > " << source.maximum() * sizeof(D) << std::endl;
		ASSERT(0);
		exit(0);
	}
	memcpy(destination, &source[0], size);
	return destination;
}

template<class D>
void* memcpy(void* destination, Array<D> source, size_t size) {
	if(size > source.maximum() * sizeof(D)) {
		std::cerr << "Size of copy exceeds source maximum: " << size
			<< " > " << source.maximum() * sizeof(D) << std::endl;
		ASSERT(0);
		exit(0);
	}
	memcpy(destination, &source[0], size);
	return destination;
}

template<class D>
void* memcpy(void* destination, ConstArray<D> source, size_t size) {
	if(size > source.maximum() * sizeof(D)) {
		std::cerr << "Size of copy exceeds source maximum: " << size
			<< " > " << source.maximum() * sizeof(D) << std::endl;
		ASSERT(0);
		exit(0);
	}
	memcpy(destination, &source[0], size);
	return destination;
}

template<class C>
Array<C> memset(Array<C> destination, int value, size_t size) {
	if(size > destination.maximum() * sizeof(C)) {
		std::cerr << "Size of set exceeds destination maximum: " << size
			<< " > " << destination.maximum() * sizeof(C) << std::endl;
		ASSERT(0);
		exit(0);
	}
	memset(&destination[0], value, size);
	return destination;
}

template<class C>
size_t fread(Array<C> destination, size_t eSize, size_t count, FILE* fp) {
	if(count * eSize > destination.maximum() * sizeof(C)) {
		std::cerr << "Size of read exceeds destination maximum: " << count * eSize
			<< " > " << destination.maximum() * sizeof(C) << std::endl;
		ASSERT(0);
		exit(0);
	}
	return fread(&destination[0], eSize, count, fp);
}

template<class C>
size_t fwrite(Array<C> source, size_t eSize, size_t count, FILE* fp) {
	if(count * eSize > source.maximum() * sizeof(C)) {
		std::cerr << "Size of write exceeds source maximum: " << count * eSize
			<< " > " << source.maximum() * sizeof(C) << std::endl;
		ASSERT(0);
		exit(0);
	}
	return fwrite(&source[0], eSize, count, fp);
}

template<class C>
size_t fwrite(ConstArray<C> source, size_t eSize, size_t count, FILE* fp) {
	if(count * eSize > source.maximum() * sizeof(C)) {
		std::cerr << "Size of write exceeds source maximum: " << count * eSize
			<< " > " << source.maximum() * sizeof(C) << std::endl;
		ASSERT(0);
		exit(0);
	}
	return fwrite(&source[0], eSize, count, fp);
}

template<class C>
void qsort(Array<C> base, size_t numElements, size_t elementSize,
		int (*compareFunction)(void const*, void const*)) {
	if(sizeof(C) != elementSize) {
		std::cerr << "Element sizes differ: " << sizeof(C) << " != " << elementSize
			<< std::endl;
		ASSERT(0);
		exit(0);
	}
	if(base.minimum() > 0 || base.maximum() < numElements) {
		std::cerr << "Array access out of bounds: " << base.minimum() << " <= 0 <= "
			<< base.maximum() << " <= " << numElements << std::endl;
		ASSERT(0);
		exit(0);
	}
	qsort(base.pointer(), numElements, elementSize, compareFunction);
}
