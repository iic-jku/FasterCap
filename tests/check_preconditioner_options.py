"""Default preconditioning must not become an unsupported explicit combination."""
from pathlib import Path
import subprocess
import sys
import tempfile

from check_capacitance import parse_last_matrix
from geometries import GEOMETRIES


def run(executable, geometry, options, expected=0):
    result = subprocess.run([executable, '-b', geometry, *options],
                            capture_output=True, text=True, timeout=10)
    output = result.stdout + result.stderr
    if result.returncode != expected:
        raise AssertionError('%r: exit %d, expected %d\n%s' %
                             (options, result.returncode, expected, output))
    if expected == 0:
        if 'Error:' in output or 'Warning:' in output:
            raise AssertionError('%r produced diagnostics\n%s' % (options, output))
        matrix = parse_last_matrix(output)
        if not matrix[1]:
            raise AssertionError('missing matrix\n' + output)
        return matrix
    if 'Capacitance matrix is:' in output or 'cannot be combined' not in output:
        raise AssertionError('unsupported combination was not rejected before solving\n' + output)


if __name__ == '__main__':
    executable = str(Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory(prefix='fc-precond-') as directory:
        geometry = str(Path(GEOMETRIES['three_cubes'].write(directory)).resolve())
        for options, reference in (([], ['-pj']), (['-pb2'], ['-pn', '-pb2']),
                                   (['-ps16'], ['-pn', '-ps16']),
                                   (['-pb2', '-pj', '-pn', '-pb2'], ['-pb2'])):
            if run(executable, geometry, options) != run(executable, geometry, reference):
                raise AssertionError('%r differs from %r' % (options, reference))
        for options in (['-pb2', '-pj'], ['-pj', '-pb2'], ['-pb2', '-ps16'],
                        ['-ps16', '-pb2'], ['-pb2', '-ps16', '-pj'], ['-pb2', '-pj', '-ap']):
            run(executable, geometry, options, 64)
