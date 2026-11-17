"""Run the ProteaseScope test suite inside Pyodide (skips slow tests)."""
import sys

import pytest

sys.path.insert(0, "/proj/src")
code = pytest.main(["/proj/tests", "-q", "-m", "not slow", "-p", "no:cacheprovider", "-x", "-s",
                    "--tb=short", "-W", "ignore::DeprecationWarning"])
print("PYTEST EXIT", int(code))
