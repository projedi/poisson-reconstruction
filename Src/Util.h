#pragma once

#include <vector>

#ifndef CPP11

#define static_assert(expr, description)
#define nullptr NULL
#define constexpr
#define override

#endif

// v.shrink_to_fit() in c++11
template<class T>
void shrink_to_fit(std::vector<T>& v) {
	std::vector<T>(v).swap(v);
}
