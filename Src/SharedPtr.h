#pragma once

#include <tr1/memory>

template<class T>
class SharedPtr {
public:
	SharedPtr(T* t): ptr_(t) { }
	T& operator*() const { return *ptr_; }
private:
	std::tr1::shared_ptr<T> ptr_;
};
