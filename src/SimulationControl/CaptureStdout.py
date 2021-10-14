try:
    from cStringIO import StringIO       # Python 2
except ImportError:
    from io import StringIO              # Python 3
import sys

# A class to capture sys.stdout
class CaptureStdout:
    def __enter__(self):
        self._stdout = sys.stdout
        self._stringio = StringIO()
        sys.stdout = self._stringio
        return self
    def __exit__(self, *args):
        self.result = self._stringio.getvalue()
        self._stringio.close()
        sys.stdout = self._stdout

# Capture the help output for the given class
def helpString(cls):
    with CaptureStdout() as output:
        help(cls)
    return output.result
