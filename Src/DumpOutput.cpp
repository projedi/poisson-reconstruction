#include <cstdarg>
#include <cstdio>
#include <cstring>

#include "DumpOutput.h"

DumpOutput& DumpOutput::instance() {
	static DumpOutput v;
	return v;
}

void DumpOutput::operator()(char const* format, ...) {
	va_list args;
	if(!noComments_) {
		va_start(args, format);
		char str[1024];
		vsprintf(str, format, args);
		if(str[strlen(str) - 1] == '\n') str[strlen(str) - 1] = 0;
		strings_.push_back(std::string(str));
		va_end(args);
	}
	if(echoStdout_) {
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
}
