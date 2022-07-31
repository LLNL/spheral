#!/usr/bin/env python

import sys
execfile("test-logs/atsr.py")
failed_tests = [t for t in state.testlist if t.status == FAILED]
if len(failed_tests) > 0:
  print("ATS failed {0} tests.".format(len(failed_tests)))
  for t in failed_tests:
    print t
  sys.exit(1)
else:
  print("ATS passed all tests.")
sys.exit(0)
