#pragma once

#include <vector>

// v.shrink_to_fit() in c++11
template<class T>
void shrink_to_fit(std::vector<T>& v) {
	std::vector<T>(v).swap(v);
}
