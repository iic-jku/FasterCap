"""Block preconditioning must preserve real, lossy and 2D physical solutions."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True

import check_capacitance
import geometries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastercap', required=True)
    args = parser.parse_args()
    executable = str(Path(args.fastercap).resolve())
    cases = [('coated_spheres_flux', ['-g'], 0.02, 0.005),
             ('coated_spheres_flux_lossy', ['-g'], 0.02, 0.005),
             ('coated_coax', [], 0.01, 0.01)]
    with tempfile.TemporaryDirectory(prefix='fc-block-') as temporary:
        root = Path(temporary)
        for name, scheme, tolerance, reciprocity in cases:
            directory = root / name
            directory.mkdir()
            geometry = geometries.GEOMETRIES[name]
            model = (geometries.write_coated_coax(str(directory), segments=512)
                     if name == 'coated_coax' else geometry.write(str(directory)))
            reference = None
            for preconditioner in ['-pj', '-pb2', '-pb32']:
                report = directory / (preconditioner + '.json')
                process = subprocess.run([executable, '-b', model, *scheme, '-m1e32',
                                          '-d1e-40', '-t1e-9', preconditioner, '-f0',
                                          '--report=' + str(report)], capture_output=True,
                                         text=True, timeout=20,
                                         env=dict(os.environ, OMP_NUM_THREADS='1'))
                if process.returncode:
                    raise AssertionError(process.stdout[-4000:] + process.stderr)
                iteration = json.loads(report.read_text())['iterations'][-1]
                assert iteration['finite'] and iteration['panels'] <= (1536 if name == 'coated_coax' else 960), (name, preconditioner, iteration)
                assert iteration['solve_status']['exit_code'] == 0
                matrix = iteration['matrix']
                rows = [[complex(a,b) for a,b in zip(real, imaginary)]
                        for real,imaginary in zip(matrix['real'], matrix['imaginary'])]
                physical = rows
                if matrix['reference_column'] is not None:
                    # Reports retain the grounded 2D reference column; the
                    # shared checker expects the reduced square matrix.
                    assert matrix['reference_column'] == geometry.conductors + 1, matrix
                    physical = [row[:geometry.conductors] for row in rows]
                errors = check_capacitance.check_matrix(physical, geometry, tolerance,
                                                        reciprocity_tolerance=reciprocity)
                if errors:
                    raise AssertionError('\n'.join(errors))
                if reference is None:
                    reference = rows
                else:
                    error = max(abs(a-b)/abs(b) for row,ref in zip(rows,reference)
                                for a,b in zip(row,ref))
                    print(name, preconditioner, 'relative difference', error)
                    assert error < 1e-6, 'block/Jacobi matrices disagree'
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        print('FAILED:', error, file=sys.stderr)
        sys.exit(1)
