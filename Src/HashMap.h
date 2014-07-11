#pragma once

#include <tr1/unordered_map>

template<class Key, class Value>
class HashMap {
public:
	typedef std::tr1::unordered_map<Key, Value> hashmap;
	typedef typename hashmap::iterator iterator;
	iterator begin() { return hashmap_.begin(); }
	iterator end() { return hashmap_.end(); }
	iterator find(Key const& key) { return hashmap_.find(key); }
	Value& operator[](Key const& key) { return hashmap_[key]; }
private:
	hashmap hashmap_;
};
