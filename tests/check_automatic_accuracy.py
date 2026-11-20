"""Automatic solves must meet physical checks and agree with a tightly solved mesh."""
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


def solve(executable, directory, model, label, options):
    report = directory / (label+'.json')
    process = subprocess.run([executable, '-b', str(model), '-g', '-pj', '-f0',
                              *options, '--report='+str(report)], cwd=directory,
                             capture_output=True, text=True, timeout=20,
                             env=dict(os.environ, OMP_NUM_THREADS='1'))
    if process.returncode:
        raise AssertionError(process.stdout[-5000:] + process.stderr)
    data = json.loads(report.read_text())
    iteration = data['iterations'][-1]
    assert data['status']['exit_code'] == 0
    assert all(it['solve_status']['exit_code'] == 0 for it in data['iterations'])
    assert iteration['finite'] and iteration['panels'] <= 2048, iteration
    matrix = iteration['matrix']
    rows = [[complex(a,b) for a,b in zip(real, imaginary)]
            for real,imaginary in zip(matrix['real'],matrix['imaginary'])]
    return rows, iteration


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastercap', required=True)
    args = parser.parse_args()
    executable = str(Path(args.fastercap).resolve())
    with tempfile.TemporaryDirectory(prefix='fc-auto-accuracy-') as temporary:
        for name, permittivity in [('moderate',10), ('high',100), ('lossy',100-20j)]:
            directory = Path(temporary) / name
            directory.mkdir()
            geometry = geometries.nested_dielectric_cubes_geometry(permittivity)
            model = geometry.write(str(directory))
            matrix, automatic = solve(executable,directory,model,'automatic',['-a0.1','-o'])
            errors = check_capacitance.check_matrix(matrix,geometry,0.1)
            if errors:
                raise AssertionError('\n'.join(errors))
            # The entire closed outer conductor surrounds the inner one.
            # Both inner shielding identities are zero in the physical matrix.
            shielding = max(abs(matrix[0][0]+matrix[0][1]),
                            abs(matrix[0][0]+matrix[1][0])) / abs(matrix[0][0])
            assert shielding <= 0.1, (name,'shielding',shielding)
            # Rebuild the same operator from the original planar input and the
            # reported refinement settings, changing only the linear tolerance.
            settings = automatic['settings']
            matched, same_operator = solve(executable,directory,model,'matched',
                ['-m'+repr(settings['mesh_refinement']),
                 '-d'+repr(settings['interaction_ratio']),
                 '-mc'+repr(settings['curvature']), '-t1e-9'])
            assert (same_operator['panels'],same_operator['links']) == (automatic['panels'],automatic['links'])
            algebraic_error = max(abs(a-b)/abs(b) for row,ref in zip(matrix,matched)
                                  for a,b in zip(row,ref))
            assert algebraic_error <= 0.025, (name,'algebraic error',algebraic_error)
            mesh = directory / 'cubes_ref.lst'
            reference, direct = solve(executable,directory,mesh,'tight',
                                      ['-m1e32','-d1e-40','-t1e-9'])
            assert direct['panels'] == automatic['panels']
            assert direct['links'] == direct['panels']*(direct['panels']-1)
            error = max(abs(a-b)/abs(b) for row,ref in zip(matrix,reference)
                        for a,b in zip(row,ref))
            print(name,'panels',automatic['panels'],'relative direct difference',error,
                  'shielding',shielding,'algebraic error',algebraic_error)
            assert error <= 0.1, (name,'direct comparison',error)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        print('FAILED:',error,file=sys.stderr)
        sys.exit(1)
