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

#include <array>
#include <cstring>
#include <string>
#include <vector>

class cmdLineReadable {
public:
	cmdLineReadable(std::string const& name): set_(false), name_(name) { }
	virtual ~cmdLineReadable() { }

	bool set() const { return set_; }
	char const* name() const { return name_.c_str(); }

	virtual int read(char** /* argv */, int /* argc */) { set_ = true; return 0; }
	virtual std::string toString() const { return ""; }
protected:
	bool set_;
	std::string name_;
};

template<class T>
class cmdLine: public cmdLineReadable {
public:
	cmdLine(std::string const& name, T const& t = T()):
		cmdLineReadable(name), value_(t) { }

	T const& value() const { return value_; }
	T& value() { return value_; }

	int read(char** argv, int argc) override;
	std::string toString() const override;
private:
	T value_;
};

// This reads the arguments in argc, matches them against "names" and sets
// the values of "r" appropriately. Parameters start with "--"
void cmdLineParse(int argc, char** argv, std::vector<cmdLineReadable*> const& r,
		bool dumpError = true);

#include "CmdLineParser.inl"
