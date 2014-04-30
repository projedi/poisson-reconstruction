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

#pragma once

#define ARRAY_DEBUG 0

#if ARRAY_DEBUG

#pragma message("[WARNING] Array debugging is enabled")

#include "Array.inl"

#define Pointer(...) Array<__VA_ARGS__>
#define ConstPointer(...) Array<__VA_ARGS__ const>

template<class C>
Array<C> NewPointer(size_t size, char const* name = nullptr)
	{ return Array<C>::New(size, name ); }

template<class C>
Array<C> AllocPointer(size_t size, char const* name = nullptr)
	{ return Array<C>::Alloc(size, false, name ); }

template<class C>
void FreePointer(Array<C>& a) { a.Free(); }

template<class C>
void DeletePointer(Array<C>& a) { a.Delete(); }

template<class C>
Array<C> NullPointer() { return Array<C>(); }

#else // !ARRAY_DEBUG

#define Pointer(...) __VA_ARGS__*
#define ConstPointer(...) __VA_ARGS__ const*

template<class C>
C* NewPointer(size_t size) { return new C[size]; }

template<class C>
C* AllocPointer(size_t size) { return (C*) malloc(sizeof(C) * size); }

#define FreePointer(...) { if(__VA_ARGS__) free(__VA_ARGS__), __VA_ARGS__ = nullptr; }

#define DeletePointer(...) { if(__VA_ARGS__) delete[] __VA_ARGS__, __VA_ARGS__ = nullptr; }

template<class C> C* NullPointer() { return nullptr; }

#endif // ARRAY_DEBUG
