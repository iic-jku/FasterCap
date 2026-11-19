"""Headless startup rejects missing input without crashing or starting a solve."""
import subprocess
import sys


def main():
    for arguments in ([], [""], [" \t\r\n"]):
        result = subprocess.run(
            [sys.argv[1], *arguments], stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, timeout=10,
        )
        if (result.returncode != 64 or "Copyright" not in result.stdout
                or "Usage:" not in result.stdout or "Running" in result.stdout):
            raise AssertionError(
                f"arguments={arguments!r}, exit={result.returncode}\n{result.stdout}"
            )


if __name__ == "__main__":
    main()
