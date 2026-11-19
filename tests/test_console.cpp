// Exercise the production tokenizer, including EOF without a %n assignment.
#include "FasterCapConsole.h"

#include <cstdio>
#include <cstring>
#include <vector>

class CConsoleParser : public FasterCapConsole
{
public:
	using FasterCapConsole::getSubstring;
};

static int failures = 0;

static void CheckToken(CConsoleParser &parser, const char *input,
		int expectedResult, const char *expectedToken, int expectedSkip)
{
	std::vector<char> token(strlen(input) + 1, 'x');
	int skip = 37;
	int result = parser.getSubstring(input, token.data(), &skip);
	if(result != expectedResult || strcmp(token.data(), expectedToken) != 0 || skip != expectedSkip) {
		std::printf("FAILED token: result=%d token=\"%s\" skip=%d; expected %d \"%s\" %d\n",
			result, token.data(), skip, expectedResult, expectedToken, expectedSkip);
		failures++;
	}
}

static void CheckEndOfTokens(CConsoleParser &parser, const char *input)
{
	std::vector<char> token(strlen(input) + 1, 'x');
	int skip = 37;
	int result = parser.getSubstring(input, token.data(), &skip);
	if(result != 1 || strcmp(token.data(), "shape.lst") != 0 || skip != 9) {
		std::printf("FAILED initial token before EOF\n");
		failures++;
		return;
	}
	const char *remaining = input + skip;
	result = parser.getSubstring(remaining, token.data(), &skip);
	if(result != EOF || token[0] != '\0' || skip != 0) {
		std::printf("FAILED EOF after token: result=%d skip=%d; expected EOF and zero\n", result, skip);
		failures++;
	}
}

int main()
{
	CConsoleParser parser;
	CheckToken(parser, "", EOF, "", 0);
	CheckToken(parser, " \t\r\n\f\v", EOF, "", 0);
	CheckToken(parser, "shape.lst", 1, "shape.lst", 9);
	CheckToken(parser, "  shape.lst ", 1, "shape.lst", 11);
	CheckToken(parser, "\"shape with spaces.lst\" -b", 1, "shape with spaces.lst", 23);
	CheckEndOfTokens(parser, "shape.lst");
	CheckEndOfTokens(parser, "shape.lst \t");

	// Do not advance the parser with an invalid count if its contract already failed.
	if(failures != 0)
		return 1;
	const char *emptyInputs[] = {"", " \t\r\n"};
	for(const char *input : emptyInputs) {
		CAutoRefGlobalVars globalVars;
		wxString error;
		if(!parser.ParseCmdLine(input, globalVars, error) ||
			!error.Contains("without enough arguments")) {
			std::printf("FAILED empty command line diagnostic\n");
			failures++;
		}
	}
	return failures != 0;
}
