#include "../stdafx.h"

#include "NumericalReport.h"
#include "../FasterCapGlobal.h"
#include "BuildProvenance.h"
#include "BuildConfiguration.h"
#include "Potential.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <locale>
#include <vector>
#include <omp.h>
#include <wx/version.h>
#include <wx/log.h>

namespace {

void WriteString(std::ostream &out, const std::string &value)
{
	const char digits[] = "0123456789abcdef";
	out << '"';
	for(std::string::const_iterator it = value.begin(); it != value.end(); ++it) {
		unsigned char c = static_cast<unsigned char>(*it);
		if(c == '"' || c == '\\') {
			out << '\\' << c;
		}
		else if(c < 0x20) {
			out << "\\u00" << digits[c >> 4] << digits[c & 15];
		}
		else {
			out << c;
		}
	}
	out << '"';
}

void WriteNumber(std::ostream &out, double value)
{
	if(value == 0.0) {
		out << (std::signbit(value) ? "-0.0" : "0.0");
	}
	else if(std::isfinite(value)) {
		out << value;
	}
	else {
		// JSON has no NaN or infinity; the record's finite flag identifies this case.
		out << "null";
	}
}

void WriteOptionalString(std::ostream &out, const char *value)
{
	if(value != NULL && value[0] != '\0') {
		WriteString(out, value);
	}
	else {
		out << "null";
	}
}

const char *StatusName(int status)
{
	switch(status) {
		case FC_NORMAL_END: return "FC_NORMAL_END";
		case FC_GENERIC_ERROR: return "FC_GENERIC_ERROR";
		case FC_COMMAND_LINE_ERROR: return "FC_COMMAND_LINE_ERROR";
		case FC_CANNOT_OPEN_FILE: return "FC_CANNOT_OPEN_FILE";
		case FC_OUT_OF_MEMORY: return "FC_OUT_OF_MEMORY";
		case FC_FILE_ERROR: return "FC_FILE_ERROR";
		case FC_CANNOT_GO_OOC: return "FC_CANNOT_GO_OOC";
		case FC_EXCEPTION_ERROR: return "FC_EXCEPTION_ERROR";
		case FC_USER_BREAK: return "FC_USER_BREAK";
		case FC_CONVERGENCE_FAILED: return "FC_CONVERGENCE_FAILED";
		case FC_INVALID_RESULT: return "FC_INVALID_RESULT";
		default: return "unknown";
	}
}

void WriteStatus(std::ostream &out, int status)
{
	out << "{\"exit_code\":" << status << ",\"name\":";
	WriteString(out, StatusName(status));
	out << '}';
}

void WriteNames(std::ostream &out, const StlStringList &names, unsigned long count)
{
	out << '[';
	StlStringList::const_iterator name = names.begin();
	for(unsigned long i = 0; i < count && name != names.end(); ++i, ++name) {
		if(i != 0) out << ',';
		WriteString(out, *name);
	}
	out << ']';
}

void WriteMatrix(std::ostream &out, const CLin_Matrix &matrix, double scale)
{
	out << '[';
	for(CLin_subscript i = 0; i < matrix.num_rows(); ++i) {
		if(i != 0) out << ',';
		out << '[';
		for(CLin_subscript j = 0; j < matrix.num_cols(); ++j) {
			if(j != 0) out << ',';
			WriteNumber(out, matrix[i][j] * scale);
		}
		out << ']';
	}
	out << ']';
}

void WriteQuality(std::ostream &out, const CLin_Matrix &real, const CLin_Matrix &imaginary, double scale)
{
	// In 2D, omit the reference-charge column, as in the ordinary square console matrix.
	const CLin_subscript n = real.num_rows();
	bool finite = real.num_cols() >= n && imaginary.num_rows() == n &&
		imaginary.num_cols() == real.num_cols();
	for(CLin_subscript i = 0; finite && i < n; ++i) {
		for(CLin_subscript j = 0; j < real.num_cols(); ++j) {
			if(!std::isfinite(real[i][j] * scale) || !std::isfinite(imaginary[i][j] * scale)) {
				finite = false;
				break;
			}
		}
	}
	out << ",\"finite\":" << (finite ? "true" : "false") << ",\"quality\":";
	if(!finite) {
		out << "null";
		return;
	}

	double norm = 0.0, defect = 0.0, worstRelative = -1.0, worstAbsolute = 0.0;
	CLin_subscript worstRow = 0, worstColumn = 0;
	unsigned long diagonal = 0, offdiagonal = 0, negativeRows = 0;
	std::vector<double> rowSums(n, 0.0);
	for(CLin_subscript i = 0; i < n; ++i) {
		if(real[i][i] <= 0.0) ++diagonal;
		for(CLin_subscript j = 0; j < n; ++j) {
			double magnitude = std::hypot(real[i][j], imaginary[i][j]);
			double difference = std::hypot(real[i][j] - real[j][i], imaginary[i][j] - imaginary[j][i]);
			norm = std::hypot(norm, magnitude);
			defect = std::hypot(defect, difference);
			rowSums[i] += real[i][j] * scale;
			if(i != j && real[i][j] >= 0.0) ++offdiagonal;
			if(j > i) {
				double pairScale = std::max(magnitude, std::hypot(real[j][i], imaginary[j][i]));
				double relative = pairScale == 0.0 ? 0.0 : difference / pairScale;
				if(relative > worstRelative) {
					worstRelative = relative;
					worstAbsolute = difference * scale;
					worstRow = i;
					worstColumn = j;
				}
			}
		}
		if(rowSums[i] < 0.0) ++negativeRows;
	}
	out << "{\"reciprocity_frobenius\":";
	if(norm == 0.0) out << "null";
	else WriteNumber(out, defect / norm);
	out << ",\"worst_reciprocal_pair\":";
	if(worstRelative < 0.0) out << "null";
	else {
		out << "{\"row\":" << worstRow + 1 << ",\"column\":" << worstColumn + 1
			<< ",\"relative_difference\":";
		WriteNumber(out, worstRelative);
		out << ",\"absolute_difference\":";
		WriteNumber(out, worstAbsolute);
		out << '}';
	}
	out << ",\"nonpositive_diagonal_count\":" << diagonal
		<< ",\"nonnegative_offdiagonal_count\":" << offdiagonal
		<< ",\"negative_row_sum_count\":" << negativeRows << ",\"row_sums_real\":[";
	for(CLin_subscript i = 0; i < n; ++i) {
		if(i != 0) out << ',';
		WriteNumber(out, rowSums[i]);
	}
	out << "]}";
}

void WriteProvenance(std::ostream &out)
{
	out << "{\"version\":";
	WriteString(out, FCG_HEADER_VERSION);
	out << ",\"compiler\":";
	WriteString(out, FC_BUILD_COMPILER);
	out << ",\"configuration\":";
	WriteString(out, FC_BUILD_CONFIGURATION);
	out << ",\"target_system\":";
	WriteString(out, FC_BUILD_SYSTEM);
	out << ",\"target_processor\":";
	WriteString(out, FC_BUILD_PROCESSOR);
	out << ",\"wxwidgets\":";
	WriteString(out, wxString(wxVERSION_STRING).ToStdString());
	out << ",\"blas_link_libraries\":";
	WriteString(out, FC_BUILD_BLAS);
	out << ",\"source_revisions\":{\"FasterCap\":";
	WriteOptionalString(out, FC_BUILD_FASTERCAP_REVISION);
	out << ",\"LinAlgebra\":";
	WriteOptionalString(out, FC_BUILD_LINALGEBRA_REVISION);
	out << ",\"Geometry\":";
	WriteOptionalString(out, FC_BUILD_GEOMETRY_REVISION);
	out << "},\"source_dirty\":{\"FasterCap\":" << FC_BUILD_FASTERCAP_DIRTY
		<< ",\"LinAlgebra\":" << FC_BUILD_LINALGEBRA_DIRTY
		<< ",\"Geometry\":" << FC_BUILD_GEOMETRY_DIRTY << "},\"omp_max_threads\":" << omp_get_max_threads();
	const char *variables[] = {"OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"};
	out << ",\"thread_environment\":{";
	for(unsigned int i = 0; i < sizeof(variables) / sizeof(variables[0]); ++i) {
		if(i != 0) out << ',';
		WriteString(out, variables[i]);
		out << ':';
		WriteOptionalString(out, std::getenv(variables[i]));
	}
	out << "}}";
}

} // namespace

CNumericalReport::CNumericalReport() : m_ulIterations(0)
{
}

bool CNumericalReport::Open(const CAutoRefGlobalVars &options)
{
	// The caller uses FasterCap's error path. wxFile's fallback logger can open
	// a modal dialog on Windows even when the solver is running in batch mode.
	wxLogNull noFileLog;
	m_stream.reset();
	m_buffer.reset();
	m_output.reset();
	if(m_file.IsOpened()) m_file.Close();
	m_ulIterations = 0;
	// Exclusive creation also protects referenced inputs and their hard/symbolic aliases.
	if(!m_file.Create(options.m_sReportFile, false)) return false;
	m_output.reset(new wxFileOutputStream(m_file));
	m_buffer.reset(new wxBufferedOutputStream(*m_output));
	m_stream.reset(new wxStdOutputStream(*m_buffer));
	m_stream->imbue(std::locale::classic());
	*m_stream << std::setprecision(std::numeric_limits<double>::max_digits10);
	*m_stream << "{\"schema_version\":1,\"input_file\":";
	WriteString(*m_stream, options.m_sFileIn);
	*m_stream << ",\"build\":";
	WriteProvenance(*m_stream);
	*m_stream << ",\"iterations\":[\n";
	return m_stream->good();
}

bool CNumericalReport::IsOpen() const
{
	return m_stream.get() != NULL;
}

bool CNumericalReport::AddIteration(const CLin_Matrix &real, const CLin_Matrix &imaginary,
	const StlStringList &names, const CAutoRefGlobalVars &options,
	unsigned char solverType, unsigned long panels, unsigned long links, int status)
{
	wxLogNull noFileLog;
	const double scale = solverType == SOLVERGLOBAL_2DSOLVER ? TWO_PI_TIMES_E0 : 1.0;
	if(m_ulIterations != 0) *m_stream << ",\n";
	*m_stream << "{\"iteration\":" << m_ulIterations << ",\"solve_status\":";
	WriteStatus(*m_stream, status);
	*m_stream << ",\"panels\":" << panels << ",\"links\":" << links << ",\"settings\":{\"automatic\":"
		<< (options.m_bAuto ? "true" : "false") << ",\"auto_tolerance\":";
	WriteNumber(*m_stream, options.m_dAutoMaxErr);
	*m_stream << ",\"gmres_tolerance\":";
	WriteNumber(*m_stream, options.m_dGmresTol);
	*m_stream << ",\"mesh_refinement\":";
	WriteNumber(*m_stream, options.m_dMeshEps);
	*m_stream << ",\"interaction_ratio\":";
	WriteNumber(*m_stream, options.m_dEpsRatio);
	*m_stream << ",\"curvature\":";
	WriteNumber(*m_stream, options.m_dMeshCurvCoeff);
	*m_stream << ",\"scheme\":";
	WriteString(*m_stream, options.m_cScheme == AUTOREFINE_GALERKIN ? "galerkin" : "collocation");
	*m_stream << ",\"preconditioner\":" << static_cast<unsigned int>(options.m_ucPrecondType) << "},\"matrix\":{\"units\":";
	WriteString(*m_stream, solverType == SOLVERGLOBAL_2DSOLVER ? "F/m" : "F");
	*m_stream << ",\"row_conductors\":";
	WriteNames(*m_stream, names, real.num_rows());
	*m_stream << ",\"column_conductors\":";
	WriteNames(*m_stream, names, real.num_cols());
	*m_stream << ",\"reference_column\":";
	if(solverType == SOLVERGLOBAL_2DSOLVER) *m_stream << real.num_cols();
	else *m_stream << "null";
	*m_stream << ",\"real\":";
	WriteMatrix(*m_stream, real, scale);
	*m_stream << ",\"imaginary\":";
	WriteMatrix(*m_stream, imaginary, scale);
	*m_stream << '}';
	WriteQuality(*m_stream, real, imaginary, scale);
	*m_stream << '}';
	++m_ulIterations;
	m_buffer->Sync();
	return m_stream->good() && m_buffer->IsOk() && m_output->IsOk();
}

bool CNumericalReport::Close(int status)
{
	wxLogNull noFileLog;
	*m_stream << "\n],\"status\":";
	WriteStatus(*m_stream, status);
	*m_stream << ",\"last_completed_iteration\":";
	if(m_ulIterations == 0) *m_stream << "null";
	else *m_stream << m_ulIterations - 1;
	*m_stream << "}\n";
	m_buffer->Sync();
	bool good = m_stream->good() && m_buffer->IsOk() && m_output->IsOk();
	m_stream.reset();
	m_buffer.reset();
	m_output.reset();
	return m_file.Close() && good;
}
