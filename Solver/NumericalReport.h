// Optional full-precision JSON report. Records are streamed, not retained in memory.
#ifndef FASTER_CAP_NUMERICAL_REPORT_H
#define FASTER_CAP_NUMERICAL_REPORT_H

#include <memory>
#include <wx/file.h>
#include <wx/wfstream.h>
#include <wx/stdstream.h>
#include "SolverGlobal.h"
#include "LinAlgebra/Mtx.h"

class CNumericalReport
{
public:
	CNumericalReport();
	bool Open(const CAutoRefGlobalVars &options);
	bool IsOpen() const;
	bool AddIteration(const CLin_Matrix &real, const CLin_Matrix &imaginary,
		const StlStringList &names, const CAutoRefGlobalVars &options,
		unsigned char solverType, unsigned long panels, unsigned long links, int status);
	bool Close(int status);

private:
	wxFile m_file;
	std::unique_ptr<wxFileOutputStream> m_output;
	std::unique_ptr<wxBufferedOutputStream> m_buffer;
	std::unique_ptr<wxStdOutputStream> m_stream;
	unsigned long m_ulIterations;
};

#endif
