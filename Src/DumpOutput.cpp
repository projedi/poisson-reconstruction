#include <cstdarg>
#include <cstring>

#include "DumpOutput.h"

DumpOutput& DumpOutput::instance() {
	static DumpOutput v;
	return v;
}

void DumpOutput::operator()(char const* format, ...) {
	va_list args;
	va_start(args, format);
	if(!noComments_) {
		char str[1024];
		vsprintf(str, format, args);
		if(str[strlen(str) - 1] == '\n') str[strlen(str) - 1] = 0;
		strings_.push_back(std::string(str));
	}
	if(echoStdout_)
		vprintf(format, args);
	va_end(args);
}
