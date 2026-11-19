"""Bounded name-overflow regressions: 20 includes and at most twelve triangles."""
from pathlib import Path
import subprocess
import sys
import tempfile


def run_case(executable, kind, name_length, rename=None, levels=20, expected=74):
    with tempfile.TemporaryDirectory(prefix='fc-names-') as temporary:
        root = Path(temporary)
        name = 'c' * name_length
        if kind == 'S':
            panels = ['S ' + name + ' 0 0 1 0', 'S ' + name + ' 1 0 1 1',
                      'S ' + name + ' 1 1 0 1', 'S ' + name + ' 0 1 0 0']
        elif kind == 'T':
            panels = ['T ' + name + ' 0 0 0 1 0 0 0 1 0']
        else:
            panels = ['Q ' + name + ' 0 0 0 1 0 0 1 1 0 0 1 0']
        if rename == 'new':
            panels.append('N ' + name + ' ' + 'r' * 480)
        elif rename == 'old':
            panels.append('N ' + 'r' * 480 + ' short')
        assert all(len(line) < 511 for line in panels)
        (root / 'panels.txt').write_text('0 panels\n' + '\n'.join(panels) + '\n')
        offset = '0 0' if kind == 'S' else '0 0 0'
        for level in range(levels):
            target = 'panels.txt' if level == levels - 1 else 'l%d.lst' % (level + 1)
            header = '2D\n' if kind == 'S' and level == 0 else ''
            (root / ('l%d.lst' % level)).write_text(header + 'C %s 1 %s\n' % (target, offset))
        result = subprocess.run([executable, '-b', '-oi', 'l0.lst'], cwd=root,
                                capture_output=True, text=True, timeout=5)
        output = result.stdout + result.stderr
        if result.returncode != expected:
            raise AssertionError('%s/%s/%d: exit %d, expected %d\n%s' %
                                 (kind, rename, name_length, result.returncode, expected, output))
        if expected == 74 and not ('name exceeds' in output and 'line ' in output
                                   and 'panels.txt' in output):
            raise AssertionError('missing file/line length diagnostic\n' + output)


if __name__ == '__main__':
    executable = str(Path(sys.argv[1]).resolve())
    for kind in ('T', 'Q', 'S'):
        run_case(executable, kind, 480)
        run_case(executable, kind, 5, rename='new')
        run_case(executable, kind, 5, rename='old')
        run_case(executable, kind, 5, expected=0)
