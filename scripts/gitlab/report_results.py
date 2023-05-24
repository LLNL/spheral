import sys

exec(compile(open("test-logs/atsr.py").read(), "test-logs/atsr.py", 'exec'))
failed_tests = [t for t in state['testlist'] if t['status'] in [FAILED,TIMEDOUT] ]
if len(failed_tests) > 0:
  print(("ATS failed {0} tests.".format(len(failed_tests))))
  for t in failed_tests:
    print(t['name'])
  sys.exit(1)
else:
  print("ATS passed all tests.")
sys.exit(0)
