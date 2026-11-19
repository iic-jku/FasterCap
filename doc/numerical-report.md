# Numerical reports

`FasterCap -b input.lst --report=result.json` writes a JSON report in addition to the existing
console output. `-e` continues to write the existing CSV format. The report does not change
solver settings, capacitances, convergence decisions, or the precision of either legacy format.

`--report=` requires one nonempty path and may appear once. Quote the argument in the shell if
the path contains spaces. The target must not exist. Exclusive creation rejects existing files, including
referenced panel files and symbolic or hard links, without changing their contents.
A report-open or report-write failure returns `FC_FILE_ERROR` (74). Command-line errors return
`FC_COMMAND_LINE_ERROR` (64) before starting a report. The shared command-line tokenizer retains dynamically sized argument storage and rejects
unterminated quotes. Report paths containing literal double quotes are rejected. Operating-system
path-length limits produce a file error; no artificial 4093-byte argument limit is imposed.
On Windows, long report paths use the absolute extended-path form, including UNC paths;
the same exclusive-creation checks still apply.

The report uses JSON schema version 1:

- `input_file`: the input path supplied to the solver
- `build`: compiler identity, actual build configuration, target system/processor, wxWidgets,
  linked BLAS libraries, the three repository revisions and dirty flags, OpenMP's maximum
  thread count, and selected thread-control environment variables
- `iterations`: one record per fully computed capacitance matrix, starting at index zero
- `status`: final solver return code and its `FC_*` name
- `last_completed_iteration`: the index of the last fully computed matrix, or `null` if none

Repository revisions and dirty states are refreshed when building. They are `null` when Git
metadata is unavailable, such as in a source archive. A dirty flag means the revision alone
does not identify all compiled source changes. Target architecture and linked-library paths
are build information; they do not claim to identify the runtime CPU model or BLAS version.
The thread environment records settings, not a measurement of actual BLAS worker counts.

Each iteration contains:

- `iteration` and `solve_status`: its index and linear-solve result, including GMRES failure
- `panels`, `links`, and `settings`: refinement and effective solver settings for that solve
- `matrix`: raw real and imaginary arrays, row/column conductor names, units, and an optional
  reference-column index
- `finite`: whether every exported matrix component is finite
- `quality`: numerical diagnostics, or `null` for a non-finite matrix

Finite double values use enough significant decimal digits to round trip, with a decimal point
independent of the process locale. Signed zero is retained. Non-finite values are written as
JSON `null`, never nonstandard `NaN` or `Infinity` tokens. Both real and imaginary arrays are
always present. Strings are escaped according to JSON syntax.

Matrix rows identify the excited conductors and columns identify the conductors whose charges
were calculated, in the same order as the console matrix. In 3D, both lists are identical. In
2D, the raw rectangular matrix also includes the last reference-conductor charge column,
which ordinary console/CSV output omits; `reference_column` is its one-based index. Values use
the same physical scaling as the console (`F` in 3D, `F/m` in 2D), with coordinates interpreted
in metres. `reference_column` is `null` in 3D. Verbose output does not change the JSON shape.

Quality diagnostics refer to the square non-reference principal matrix in 2D:

- `reciprocity_frobenius` is `||C - C^T||F / ||C||F`, using the transpose also for complex
  capacitances; it is `null` for a zero matrix
- `worst_reciprocal_pair` identifies the one-based row and column with maximum
  `|Cij - Cji| / max(|Cij|, |Cji|)` and gives both that ratio and the absolute difference;
  two zero entries give zero difference, and a one-conductor matrix has no pair (`null`)
- `nonpositive_diagonal_count` and `nonnegative_offdiagonal_count` count real-part sign checks
- `row_sums_real` contains sums of the real parts, and `negative_row_sum_count` counts sums
  strictly below zero

These are raw diagnostics without acceptance tolerances. A tiny coupling can have a large
relative pair discrepancy, and a slightly negative row sum can arise from numerical error.
Consumers must choose physically justified absolute/relative budgets; a low reciprocity defect
alone does not establish capacitance accuracy. No averaging or clipping is applied.

The last completed matrix can be unconverged or invalid. Always inspect `status`, its
`solve_status`, and `finite`; do not interpret the presence of a matrix as success. If input
reading or another stage fails after earlier iterations, those earlier records remain available
with the failure status. A partial excitation sweep is omitted. `-oi` performs no capacitance
solve and therefore produces an empty iteration list. An abruptly killed process or failed disk
write can leave an incomplete report; JSON parse failure must not be treated as success.

Records are streamed to disk rather than retained as a second history in memory. Every
iteration includes its matrix so a consumer can independently compare entries and convergence
between runs without relying on rounded console values.
