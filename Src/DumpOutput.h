#pragma once

#include <string>
#include <vector>

class DumpOutput {
public:
	static DumpOutput& instance();
	void setEchoStdout(bool v) { echoStdout_ = v; }
	void setNoComments(bool v) { noComments_ = v; }
	std::vector<std::string> strings() const { return strings_; }
	void operator()(char const* format, ...);
	void resetStrings(std::vector<std::string> const& strs) { strings_ = strs; }
private:
	DumpOutput(): echoStdout_(false), noComments_(false) { }
private:
	bool echoStdout_;
	bool noComments_;
	std::vector<std::string> strings_;
};
