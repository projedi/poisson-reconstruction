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

#include <iostream>
#include <vector>

#include "SharedPtr.h"

struct AllocatorState {
	AllocatorState(): index(0), remains(0) { }
	size_t index;
	size_t remains;
};

/** This templated class assists in memory allocation and is well suited for instances
  * when it is known that the sequence of memory allocations is performed in a stack-based
  * manner, so that memory allocated last is released first. It also preallocates memory
  * in chunks so that multiple requests for small chunks of memory do not require separate
  * system calls to the memory manager.
  * The allocator is templated off of the class of objects that we would like it to allocate,
  * ensuring that appropriate constructors and destructors are called as necessary.
  */
template<class T>
class Allocator {
public:
	Allocator(): block_size_(0) { }
	~Allocator() { reset(); }

	/** This method returns the memory state of the allocator. */
	AllocatorState getState() const { return state_; }

	/** This method initiallizes the constructor and the block_size variable specifies the
	  * the number of objects that should be pre-allocated at a time. */
	void set(size_t block_size);

	/** This method returns a pointer to an array of elements objects. If there is
	  * left over pre-allocated memory, this method simply returns a pointer to
	  * the next free piece of memory, otherwise it pre-allocates more memory.
	  * Note that if the number of objects requested is larger than the value
	  * block_size_ with which the allocator was initialized, the request for
	  * memory will fail.
	  */
	T* newElements(size_t elements = 1);
private:
	/** This method is the allocators destructor. It frees up any of the memory that
	  * it has allocated. */
	void reset();
private:
	size_t block_size_;
	AllocatorState state_;
	// SharedPtr so that inner vector does not copy when memory_ resizes.
	std::vector<SharedPtr<std::vector<T>>> memory_;
};

template<class T>
void Allocator<T>::reset() {
	memory_.clear();
	state_ = AllocatorState();
	block_size_ = 0;
}

template<class T>
void Allocator<T>::set(size_t block_size) {
	reset();
	block_size_ = block_size;
	state_.remains = block_size_;
	memory_.push_back(new std::vector<T>(block_size_));
}

template<class T>
T* Allocator<T>::newElements(size_t elements) {
	if(!elements) return nullptr;
	if(elements > block_size_) {
		std::cerr << "Allocator Error, elements bigger than block size:"
					 << elements << " > " << block_size_ << std::endl;
		return nullptr;
	}
	if(state_.remains < elements) {
		if(state_.index == memory_.size() - 1)
			memory_.push_back(new std::vector<T>(block_size_));
		++state_.index;
		state_.remains = block_size_;
	}
	T* mem = &((*memory_[state_.index])[block_size_ - state_.remains]);
	state_.remains -= elements;
	return mem;
}
