"""Invalid user text must be printed literally, never interpreted by printf."""
import subprocess
import sys

executable = sys.argv[1]
for option in ('-x%s%s%s', '-x%n', '-x%p', '-x%%', '-m0.1%'):
    result = subprocess.run([executable, '-b', option, 'unused.lst'],
                            capture_output=True, text=True, timeout=5)
    output = result.stdout + result.stderr
    if result.returncode != 64 or option not in output:
        raise SystemExit('invalid option was not reported literally: ' + repr(option)
                         + '\nexit=' + str(result.returncode) + '\n' + output)
