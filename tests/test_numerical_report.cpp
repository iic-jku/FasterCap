// Write synthetic reports for Python's independent JSON parser and numerical checks.
#include <cmath>
#include <cstdio>
#include <limits>
#include <locale>
#include <string>
#include <wx/log.h>

#include "NumericalReport.h"

class CCommaDecimal : public std::numpunct<char>
{
protected:
	char do_decimal_point() const { return ','; }
};

class CReportLog : public wxLog
{
public:
	unsigned long m_ulMessages = 0;
protected:
	void DoLogText(const wxString &) override { ++m_ulMessages; }
};

static bool writeReport(const std::string &directory, const char *filename,
	const CLin_Matrix &real, const CLin_Matrix &imaginary, unsigned char dimension, int status)
{
	CAutoRefGlobalVars options;
	options.m_sFileIn = "input\\with\"quotes\tand\ncontrols.lst";
	options.m_sReportFile = directory + "/" + filename;
	StlStringList names;
	names.push_back("first\\\"\t\n\r\b\f\x01");
	names.push_back("second");
	CNumericalReport report;
	return report.Open(options) && report.AddIteration(real, imaginary, names, options,
		dimension, 12, 34, status) && report.Close(status);
}

int main(int argc, char **argv)
{
	if(argc != 2) return 1;
	std::locale::global(std::locale(std::locale::classic(), new CCommaDecimal));
	CLin_Matrix real(2, 2), imaginary(2, 2, 0.0);
	real[0][0] = std::nextafter(1.0, 2.0);
	real[0][1] = -0.12345678901234567;
	real[1][0] = -0.125;
	real[1][1] = std::numeric_limits<double>::min();
	imaginary[1][1] = -0.0;
	if(!writeReport(argv[1], "real.json", real, imaginary, SOLVERGLOBAL_3DSOLVER, FC_NORMAL_END)) return 1;
	imaginary[0][0] = -std::nextafter(0.25, 1.0);
	imaginary[0][1] = 0.0625;
	imaginary[1][0] = 0.03125;
	imaginary[1][1] = -0.25;
	if(!writeReport(argv[1], "complex.json", real, imaginary, SOLVERGLOBAL_3DSOLVER, FC_CONVERGENCE_FAILED)) return 1;

	real.newsize(1, 2);
	imaginary.newsize(1, 2);
	real[0][0] = 2.0;
	real[0][1] = -2.0;
	imaginary[0][0] = -0.5;
	imaginary[0][1] = 0.5;
	if(!writeReport(argv[1], "2d.json", real, imaginary, SOLVERGLOBAL_2DSOLVER, FC_NORMAL_END)) return 1;

	real[0][0] = std::numeric_limits<double>::quiet_NaN();
	real[0][1] = std::numeric_limits<double>::infinity();
	imaginary[0][0] = -std::numeric_limits<double>::infinity();
	if(!writeReport(argv[1], "nonfinite.json", real, imaginary, SOLVERGLOBAL_2DSOLVER, FC_INVALID_RESULT)) return 1;

	CAutoRefGlobalVars options;
	options.m_sFileIn = "missing.lst";
	options.m_sReportFile = std::string(argv[1]) + "/empty.json";
	CNumericalReport report;
	if(!report.Open(options) || !report.Close(FC_CANNOT_OPEN_FILE)) return 1;
	// The caller reports file errors through FasterCap's error path. wxWidgets'
	// default fallback can display a modal dialog even in a Windows batch run.
	options.m_sReportFile = std::string(argv[1]) + "/missing/report.json";
	CReportLog log;
	wxLog *previous = wxLog::SetActiveTarget(&log);
	bool enabled = wxLog::EnableLogging(true);
	bool opened = report.Open(options);
	wxLog::FlushActive();
	bool restored = wxLog::IsEnabled();
	wxLog::EnableLogging(enabled);
	wxLog::SetActiveTarget(previous);
	if(opened || log.m_ulMessages != 0 || !restored) {
		std::fprintf(stderr, "report file error escaped to wxWidgets logging\n");
		return 1;
	}
	return 0;
}
