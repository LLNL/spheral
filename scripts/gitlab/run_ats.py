#!/usr/bin/env python3

import sys, subprocess, argparse, os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from spheralutils import sexe

# If the number of failed tests exceeds this value, ATS is not rerun
max_test_failures = 10
# Number of times to rerun the ATS tests
max_reruns = 1

#------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()

    # Spec args
    parser.add_argument('--test-alloc', type=str, nargs="+",
                        help='Allocation command for the machine.')
    parser.add_argument('--ats-file', type=str,
                        help='ATS test file to run.')
    parser.add_argument('--ci-build-dir', type=str,
                        help='CI build directory.')
    parser.add_argument('--ci-install-dir', type=str,
                        default="build_gitlab/install",
                        help="Location of Spheral installation "+\
                        "relative to --ci-build-dir")
    return parser.parse_args()

#------------------------------------------------------------------------------

# Run ats.py to check results and return the number of failed tests
def report_results(output_dir):
    ats_py = os.path.join(output_dir, "atsr.py")
    if (not os.path.exists(ats_py)):
        print(f"{ats_py} does not exists")
        sys.exit(1)
    exec(compile(open(ats_py).read(), ats_py, 'exec'), globals())
    state = globals()["state"]
    failed_tests = [t for t in state['testlist'] if t['status'] in [FAILED,TIMEDOUT] ]
    if len(failed_tests) > 0:
        print(f"ATS failed {len(failed_tests)} tests.")
        for t in failed_tests:
            print(t['name'])
        return len(failed_tests)
    else:
        print("ATS passed all tests.")
        return 0

#------------------------------------------------------------------------------

# Run the tests and check if any failed
def run_and_report(run_command, ci_output, num_runs):
    if (num_runs > max_reruns):
        print("Exceeded number of ATS reruns")
        sys.exit(1)
    sexe(run_command)
    tests_passed = report_results(ci_output)
    if (tests_passed == 0):
        if (num_runs > 0):
            print("WARNING: Some tests were run multiple times")
        sys.exit(0)
        # This should be added back in once Jacamar can handle exit codes properly
        # if (num_runs == 0):
        #     sys.exit(0)
        # else:
        #     sys.exit(80)
    elif (tests_passed >= max_test_failures):
        print("Too many test failures, not rerunning ATS")
        sys.exit(1)
    else:
        rerun_command = run_command
        if (num_runs == 0):
            ats_cont_file = os.path.join(ci_output, "continue.ats")
            if (not os.path.exists(ats_cont_file)):
                print(f"{ats_cont_file} not found, ATS cannot be rerun")
                sys.exit(1)
            rerun_command = f"{run_command} {ats_cont_file}"
        print("WARNING: Test failure, rerunning ATS")
        run_and_report(rerun_command, ci_output, num_runs + 1)

#------------------------------------------------------------------------------

def run_ats_test(args):
    build_gl_dir = os.path.join(args.ci_build_dir, args.ci_install_dir)
    ats_file = os.path.join(build_gl_dir, args.ats_file)
    if (not os.path.exists(ats_file)):
        print(f"{ats_file} does not exists")
        sys.exit(1)
    lcats_test = os.path.join(build_gl_dir, "spheral-lcatstest")
    if (not os.path.exists(lcats_test)):
        print(f"{lcats_test} does not exists")
    ats_configs = ' --timelimit="45m"'
    test_alloc = " ".join(args.test_alloc)
    run_command = f"{test_alloc} {lcats_test} --logs test-logs {ats_file} {ats_configs}"
    ci_output = os.path.join(args.ci_build_dir, "test-logs")
    run_and_report(run_command, ci_output, 0)

#------------------------------------------------------------------------------

def main():
    args = parse_args()
    run_ats_test(args)

if __name__ == "__main__":
  main()
