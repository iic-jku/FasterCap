/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


#include "FasterCapConsole.h"

#include "Solver/SolveCapacitance.h"

#include <cmath>
#include <csignal>
#include <vector>

// Ctrl-C or kill: stop the run as the stop button of the GUI does, so that it deletes its
// out-of-core temporary files; a second signal terminates at once
static void StopOnSignal(int sig)
{
	g_bFCContinue = false;
	signal(sig, SIG_DFL);
}

#ifdef FCG_HEADLESS

#include <iostream>
#include <wx/debug.h>

// failed asserts are printed, not shown in a dialog: wxWidgets opens one on Windows,
// also in console programs, and nobody closes it in a batch run
static void ConsoleAssertHandler(const wxString& file, int line, const wxString& func,
                                 const wxString& cond, const wxString& msg)
{
    std::cerr << file << "(" << line << "): assert \"" << cond << "\" failed in "
              << func << "(): " << msg << std::endl;
}

int main(int argc, char **argv)
{
    int i, ret;
    char bOption;

    wxSetAssertHandler(ConsoleAssertHandler);

    // the character after '-b' ('-b?', '-bv'), as FasterCapApp::OnInit() finds it for the GUI
    bOption = '\0';
    for(i=1; i<argc; i++) {
        if(argv[i][0] == '-' && argv[i][1] == 'b') {
            bOption = argv[i][2];
        }
    }

    FasterCapConsole console;

    ret = console.main(argc, argv, bOption);

    return ret;
}

int FasterCapConsole::main(int& argc, char **argv, char bOption)
#else
int FasterCapConsole::main(int& argc, wxChar **argv, char bOption)
#endif // FCG_HEADLESS
{
    bool ret;
	CSolveCap solveMain;
	int i, retStatus;
	wxString inputArgs, errMsg;
	CAutoRefGlobalVars globalVars, defGlobalVars;

	// copy input arguments in a string for proper parsing (e.g. taking care of file names with spaces)
	inputArgs.Clear();
	for(i=1; i<argc; i++) {
		inputArgs += argv[i];
		// this puts an extra trailing space in the end but that's ok anyway for the parser
		inputArgs += wxT(" ");
	}

	// dump command args for debug
	//printf("%s\n", (const char*)inputArgs);

	// parse command line
	ret = ParseCmdLine((const char*)inputArgs, globalVars, errMsg);

	// if no parser error
	if(ret == false && bOption != '?' && bOption != 'v') {
		// actually run FasterCap
		signal(SIGINT, StopOnSignal);
		signal(SIGTERM, StopOnSignal);
		retStatus = FC_NORMAL_END;
		retStatus = solveMain.Run(globalVars);
	}
	else {
		if(bOption == '?' || bOption == 'v') {
			retStatus = FC_NORMAL_END;
			errMsg = wxT("");
		}
		else {
			retStatus = FC_COMMAND_LINE_ERROR;
		}

		LogMsg(FCG_HEADER_VERSION);
		LogMsg("\n");
		LogMsg(FCG_HEADER_COPYRIGHT);
		LogMsg(" ");
		LogMsg(FCG_HEADER_WEBSITE);
		LogMsg("\n");

		if(bOption != 'v') {
			LogMsg(FCG_HEADER_FORK);
			LogMsg("\n\n");

			// print error
			ErrMsg((const char*)errMsg);
			LogMsg("Usage: %s <input file> [-a<relative error>] [-ap]\n", (const char*)argv[0]);
			LogMsg("                 [-m<mesh>] [-mc<mesh curvature] [-t<tolerance>]\n");
			LogMsg("                 [-d<interaction coeff>] [-f<outofcore>] [-g]\n");
			LogMsg("                 [-pn] [-pj] [-ps<dimension>]\n");
			LogMsg("                 [-o] [-oi] [-e] [-r] [-c] [-i] [-v]\n");
			LogMsg("                 [-b|-b?|-bv]\n");
			LogMsg("DEFAULT VALUES:\n");
			LogMsg("  -a:  Automatically calculate settings, stop when\n");
			LogMsg("       relative error is lower than <relative error>, e.g. 0.01\n");
			LogMsg("  -ap: Automatic preconditioner usage\n");
			LogMsg("  -m:  Mesh relative refinement value = %g\n", defGlobalVars.m_dMeshEps);
			LogMsg("  -mc: Mesh curvature coefficient = %g\n", defGlobalVars.m_dMeshCurvCoeff);
			LogMsg("  -t:  GMRES iteration tolerance = %g\n", defGlobalVars.m_dGmresTol);
			LogMsg("  -d:  Direct potential interaction coefficient to mesh refinement ratio = %g\n", defGlobalVars.m_dEpsRatio);
			LogMsg("  -f:  Out-Of-Core free memory to link memory condition = %g\n", defGlobalVars.m_dOutOfCoreRatio);
			LogMsg("  -g:  Use Galerkin scheme\n");
			LogMsg("  -pn: Use no preconditioner\n");
			LogMsg("  -pj: Use Jacobi Preconditioner\n");
			LogMsg("  -ps: Use two-levels preconditioner with dimension = %d\n", defGlobalVars.m_uiSuperPreDim);
			LogMsg("OPTIONS:\n");
			LogMsg("  -o:  Output refined geometry in FastCap2 format\n");
			LogMsg("  -oi: Dump input geometry in FasterCap format and stop\n");
			LogMsg("  -e:  Output capacitance matrix to file\n");
			LogMsg("  -r:  Dump Gmres residual at each iteration\n");
			LogMsg("  -c:  Dump charge densities in output file\n");
			LogMsg("  -i:  Dump detailed time and memory information\n");
			LogMsg("  -v:  Verbose output\n");
			LogMsg("  -b:  Launch as console/shell application without GUI\n");
			LogMsg("  -b?: Print console usage (this text)\n");
			LogMsg("  -bv: Print only the version\n");
#ifndef __WXMSW__
			LogMsg("  -h: Open only the help system\n");
#endif
		}
	}

	// if return code is an error, print it
	solveMain.PrintRetError(retStatus);

	return retStatus;
}


// reads the value of a numeric option from 'valueStr' into 'value';
// returns 'true' if it is no finite number from 'min' (excluded unless 'minIncluded') to 'max' (excluded),
// or if other characters follow the number
static bool BadNumber(const char *valueStr, double *value, double min, bool minIncluded, double max)
{
	int length;

	if(sscanf(valueStr, "%lf%n", value, &length) != 1 || valueStr[length] != '\0' || std::isfinite(*value) == false) {
		return true;
	}
	return *value < min || (*value == min && minIncluded == false) || *value >= max;
}

// reads the value of an integer option from 'valueStr' into 'value';
// returns 'true' if it is no whole number from 'min' to 'max' (both included), as BadNumber() checks it
static bool BadInteger(const char *valueStr, unsigned int *value, unsigned int min, unsigned int max)
{
	double number;

	if(BadNumber(valueStr, &number, min, true, max + 1.0) == true || number != std::floor(number)) {
		return true;
	}
	*value = (unsigned int)number;
	return false;
}

bool FasterCapConsole::ParseOptionValue(const char *option, const char *valueStr, CAutoRefGlobalVars &globalVars, wxString &errMsg)
{
	wxString what, range;
	bool bad;

	if(strcmp(option, "a") == 0) {
		bad = BadNumber(valueStr, &(globalVars.m_dAutoMaxErr), 0.0, false, 1.0);
		what = wxT("max error value for auto option");
		range = wxT("a number > 0 and < 1");
	}
	else if(strcmp(option, "t") == 0) {
		bad = BadNumber(valueStr, &(globalVars.m_dGmresTol), 0.0, false, 1.0);
		what = wxT("GMRES iteration tolerance");
		range = wxT("a number > 0 and < 1");
	}
	else if(strcmp(option, "m") == 0) {
		bad = BadNumber(valueStr, &(globalVars.m_dMeshEps), 0.0, false, HUGE_VAL);
		what = wxT("mesh relative refinement value");
		range = wxT("a number > 0");
	}
	else if(strcmp(option, "mc") == 0) {
		bad = BadNumber(valueStr, &(globalVars.m_dMeshCurvCoeff), 1.0, true, HUGE_VAL);
		what = wxT("mesh curvature coefficient");
		range = wxT("a number >= 1");
	}
	else if(strcmp(option, "d") == 0) {
		bad = BadNumber(valueStr, &(globalVars.m_dEpsRatio), 0.0, false, HUGE_VAL);
		what = wxT("direct potential interaction coefficient to mesh refinement ratio value");
		range = wxT("a number > 0");
	}
	else if(strcmp(option, "f") == 0) {
		bad = BadNumber(valueStr, &(globalVars.m_dOutOfCoreRatio), 0.0, true, HUGE_VAL);
		what = wxT("Out-Of-Core free memory to link memory condition");
		range = wxT("a number >= 0");
	}
	else if(strcmp(option, "pb") == 0) {
		bad = BadInteger(valueStr, &(globalVars.m_uiBlockPreSize), SOLVE_MIN_BLOCK_PRECOND_NUM, SOLVE_MAX_BLOCK_PRECOND_NUM);
		what = wxT("preconditioner dimension");
		range = wxString::Format(wxT("a whole number from %d to %d"), SOLVE_MIN_BLOCK_PRECOND_NUM, SOLVE_MAX_BLOCK_PRECOND_NUM);
	}
	else if(strcmp(option, "ps") == 0) {
		bad = BadInteger(valueStr, &(globalVars.m_uiSuperPreDim), SOLVE_MIN_SUPER_PRECOND_NUM, SOLVE_MAX_SUPER_PRECOND_NUM);
		what = wxT("preconditioner dimension");
		range = wxString::Format(wxT("a whole number from %d to %d"), SOLVE_MIN_SUPER_PRECOND_NUM, SOLVE_MAX_SUPER_PRECOND_NUM);
	}
	else {
		errMsg = wxString::Format(wxT("unknown numeric option '-%s'"), option);
		return true;
	}

	if(bad == true) {
		errMsg = wxString::Format(wxT("bad %s '%s', must be %s"), what, valueStr, range);
	}
	return bad;
}

// returns 'true' if errors in the parsing
bool FasterCapConsole::ParseCmdLine(const char *commandStr, CAutoRefGlobalVars &globalVars, wxString &errMsg)
{
	// no argument is longer than the command line
	std::vector<char> argBuffer(strlen(commandStr) + 1);
	char *argStr = argBuffer.data(), *cmdStr;
	int res, skip;
	bool cmderr;

    // clean error message
    errMsg.Clear();

	// copy pointer, not to change 'commandStr'
	// (required when the function is called by automation passing a 'const')
	cmdStr = (char *)commandStr;

	// defaults

	globalVars.Reset();

	// parse input arguments
	//

	res = getSubstring(cmdStr, argStr, &skip);
	cmdStr += skip;
	if( res == 0 || res == EOF)  {
		errMsg = wxString::Format(wxT("FasterCap launched without enough arguments! You must specify at least the input file name\n"));
		errMsg += wxString::Format(wxT("Command line is: %s\n"), cmdStr);
		return true;
	}

	cmderr = false;
	while(res != 0 && res != EOF && cmderr != true) {
		if(argStr[0] == '-') {

			// '-m' is mesh relative refinement level
			// '-mc' is the curvature coefficient for gradient meshing near sharp edges
			if(argStr[1] == 'm') {
				if(argStr[2] == 'c') {
					cmderr = ParseOptionValue("mc", &(argStr[3]), globalVars, errMsg);
				}
				else {
					cmderr = ParseOptionValue("m", &(argStr[2]), globalVars, errMsg);
				}
			}

			// '-d' is direct potential interaction coefficient to mesh refinement ratio
			else if(argStr[1] == 'd') {
				cmderr = ParseOptionValue("d", &(argStr[2]), globalVars, errMsg);
			}

			// '-f' is Out-Of-Core free memory to link memory condition (0 = don't go OOC)
			else if(argStr[1] == 'f') {
				cmderr = ParseOptionValue("f", &(argStr[2]), globalVars, errMsg);
			}

			// '-p' is type of preconditioner
			else if(argStr[1] == 'p') {
				if(argStr[2] == 'N' || argStr[2] == 'n') {
					globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_NONE;
				}
				else if (argStr[2] == 'J' || argStr[2] == 'j') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_JACOBI;
				}
				else if (argStr[2] == 'B' || argStr[2] == 'b') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_BLOCK;
					cmderr = ParseOptionValue("pb", &(argStr[3]), globalVars, errMsg);
				}
				else if (argStr[2] == 'S' || argStr[2] == 's') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_SUPER;
					cmderr = ParseOptionValue("ps", &(argStr[3]), globalVars, errMsg);
				}
				else {
					cmderr = true;
					errMsg = wxString::Format(wxT("bad preconditioner type '%s'"), &argStr[2]);
				}
			}

			// '-c' is dump charge densities in output file
			else if(argStr[1] == 'c') {
                globalVars.m_bOutputCharge = true;
			}
			
			// '-r' is dump residual at every GMRES step
			else if(argStr[1] == 'r') {
				globalVars.m_bDumpResidual = true;
			}

			// '-i' is dump detailed time and memory information
			else if(argStr[1] == 'i') {
				globalVars.m_bDumpTimeMem = true;
			}

			// '-v' is provide verbose output
			else if(argStr[1] == 'v') {
				globalVars.m_bVerboseOutput = true;
			}

			// '-t' is GMRES tolerance to stop iteration
			else if(argStr[1] == 't') {
				cmderr = ParseOptionValue("t", &(argStr[2]), globalVars, errMsg);
			}

			// '-o' is output refined geometry to file (in FastCap2 format),
			// '-oi' is dump input geometry, no refinement, in FastCap2 format and stop
			//       (so actual value of 'm_bOutputGeo' does not matter)
			else if(argStr[1] == 'o') {
				if(argStr[2] == 'i') {
					globalVars.m_bDumpInputGeo = true;
				}
				else {
                    globalVars.m_bOutputGeo = true;
				}
			}

			// '-e' is output capacitance matrix to file
			else if(argStr[1] == 'e') {
				globalVars.m_bOutputCapMtx = true;
			}

			// '-g' is galerkin scheme (as opposite to default collocation)
			else if(argStr[1] == 'g') {
				globalVars.m_cScheme = AUTOREFINE_GALERKIN;
			}

			// '-a' is automatic refinement and selection of parameters
			else if(argStr[1] == 'a') {
				if(argStr[2] == 'p' || argStr[2] == 'P') {
					globalVars.m_bAutoPrecond = true;
				}
				else {
					globalVars.m_bAuto = true;
					cmderr = ParseOptionValue("a", &(argStr[2]), globalVars, errMsg);
				}
			}

			// '-b' is the console option; not to be passed to FasterCap engine
			else if(argStr[1] != 'b') {
				errMsg = wxString::Format(wxT("illegal option -- '-%s'"), &(argStr[1]));
				cmderr = true;
			}
		}
		// isn't an option, must be the input file
		else {
			globalVars.m_sFileIn = argStr;
		}

		// get next argument
		res = getSubstring(cmdStr, argStr, &skip);
		cmdStr += skip;
	}

	if(cmderr == true) {
		errMsg = wxString::Format(wxT("%s: %s\n"), commandStr, errMsg);
	}

	return cmderr;
}

// Special version of sscanf which handles also file names with spaces inside,
// provided they are surrounded by '"'.
// Note that this routine is specialized to retrieve only strings;
// 'substr' must hold strlen(buffer) + 1 characters
int FasterCapConsole::getSubstring(const char *buffer, char *substr, int *skip)
{
	int res, openPos, deltaClosePos, startPos;
	std::vector<char> tmpBuffer(strlen(buffer) + 1);
	char *openPosPtr, *closePosPtr, *tmpStr = tmpBuffer.data();

	substr[0] = '\0';

	// read a piece of string
	res = sscanf(buffer, "%s%n", tmpStr, skip);
	// if not finished
	if( res != EOF ) {
		// find if it does contain any '"'
		openPosPtr = (char *)strchr(buffer, '"');
        // if found and if within tmpStr (before any space)
		openPos = (int)(openPosPtr - buffer);
		if( openPosPtr != NULL && openPos < (int)strlen(tmpStr) ) {
			// search for the closing '"'
			closePosPtr = (char *)strchr(buffer + openPos + 1, '"');
			// if no closing '"', assume the end is the end of the
			// string in the buffer, regardless of any spaces or tabs
			if(closePosPtr == NULL) {
				// build 'substr' string by collating the different pieces,
				// skipping the '"'
				strncpy(tmpStr, buffer, openPos);
				tmpStr[openPos] = '\0';
				strcat(tmpStr, buffer + openPos + 1);
				// remove any trailing space left
				startPos = strspn(tmpStr, " \t\n");
				strcpy(substr, tmpStr + startPos);
				*skip = strlen(buffer);
			}
			else {
				// build 'substr' string by collating the different pieces,
				// skipping the '"'s
				deltaClosePos = (int)(closePosPtr - openPosPtr);
				strncpy(tmpStr, buffer, openPos);
				tmpStr[openPos] = '\0';
				strncat(tmpStr, buffer + openPos + 1, deltaClosePos - 1);
				// remove any trailing space left
				startPos = strspn(tmpStr, " \t\n");
				strcpy(substr, tmpStr + startPos);
				*skip = openPos + deltaClosePos + 1;
			}
		}
		// if no '"' at all, simply copy scanned string
		else {
			strcpy(substr, tmpStr);
		}
	}

	return res;
}

