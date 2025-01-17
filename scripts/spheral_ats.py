#!/usr/bin/env python3

import os, time, sys, subprocess, argparse
import ats.util.generic_utils as ats_utils
import SpheralConfigs
import mpi

# This is a wrapper for running Spheral through ATS

# Options for running CI
# If the number of failed tests exceeds this value, ATS is not rerun
max_test_failures = 10
# Number of times to rerun the ATS tests
max_reruns = 1

cur_dir = os.path.dirname(__file__)
# Set current directory to install prefix
if (os.path.islink(__file__)):
    cur_dir = os.path.join(cur_dir, os.readlink(__file__))
install_prefix = os.path.join(cur_dir, "..")
ats_exe = os.path.join(install_prefix, ".venv/bin/ats")
spheral_exe = os.path.join(install_prefix, "spheral")

# Benchmark file directory
# This is passed into both ATS and Caliper
benchmark_dir = "/usr/WS2/sduser/Spheral/benchmarks"

#------------------------------------------------------------------------------
# Run ats.py to check results and return the number of failed tests
def report_results(output_dir):
    ats_py = os.path.join(output_dir, "atsr.py")
    if (not os.path.exists(ats_py)):
        raise Exception("ats.py does not exists. Tests likely did not run.")
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
        raise Exception ("Exceeded number of ATS reruns")
    ats_cont_file = os.path.join(ci_output, "continue.ats")
    new_run_command = run_command
    if (os.path.exists(ats_cont_file) and num_runs == 0):
        new_run_command = f"{run_command} {ats_cont_file}"
        print("Restarting from previous job")
    try:
        subprocess.run(new_run_command, shell=True, check=True, text=True)
    except Exception as e:
        print(e)
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
        raise Exception("Too many test failures, not rerunning ATS")
    else:
        rerun_command = run_command
        if (num_runs == 0):
            ats_cont_file = os.path.join(ci_output, "continue.ats")
            if (not os.path.exists(ats_cont_file)):
                raise Exception(f"{ats_cont_file} not found, ATS cannot be rerun")
            rerun_command = f"{run_command} {ats_cont_file}"
        print("WARNING: Test failure, rerunning ATS")
        run_and_report(rerun_command, ci_output, num_runs + 1)

#------------------------------------------------------------------------------
# Add any build specific ATS arguments
def install_ats_args():
    install_args = []
    if (SpheralConfigs.build_type() == "Debug"):
        install_args.append("--level 99")
    if (mpi.is_fake_mpi()):
        install_args.append("--filter='np<2'")
    comp_configs = SpheralConfigs.component_configs()
    test_comps = ["FSISPH", "GSPH", "SVPH"]
    for ts in test_comps:
        if ts not in comp_configs:
            install_args.append(f"--filter='not {ts.lower()}'")
    return install_args

#---------------------------------------------------------------------------
# Main routine
#---------------------------------------------------------------------------
def main():
    test_log_name = "test-logs"
    toss_machine_names = ["rzgenie", "rzwhippet", "rzhound", "ruby"]
    blueos_machine_names = ["rzansel", "lassen"]
    ci_launch_flags = {"ruby": "--res=ci", "lassen": "-q pci"}
    temp_uname = os.uname()
    hostname = temp_uname[1]
    sys_type = os.getenv("SYS_TYPE")
    # Use ATS to for some machine specific functions
    if "MACHINE_TYPE" not in os.environ:
        ats_utils.set_machine_type_based_on_sys_type()

    #---------------------------------------------------------------------------
    # Setup argument parser
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     usage="""
                                     ./spheral-ats --numNodes 2 tests/integration.ats --filter="level<100"
                                     """,
                                     description="""
                                     Launches and runs Spheral using the ATS system.
                                     Must provide an ATS file (either python or .ats).
                                     Any unrecognized arguments are passed as inputs to the ATS file.
                                     """)
    parser.add_argument("--numNodes", type=int,
                        default=None,
                        help="Number of nodes to allocate.")
    parser.add_argument("--timeLimit", type=int,
                        default=None,
                        help="Time limit for allocation.")
    parser.add_argument("--ciRun", action="store_true",
                        help="Option to only be used by the CI")
    parser.add_argument("--atsHelp", action="store_true",
                        help="Print the help output for ATS. Useful for seeing ATS options.")
    parser.add_argument("--threads", type=int, default=None,
                        help="Set number of threads per rank to use. Only used by performance.py")
    options, unknown_options = parser.parse_known_args()
    if (options.atsHelp):
        subprocess.run(f"{ats_exe} --help", shell=True, check=True, text=True)
        return

    #---------------------------------------------------------------------------
    # Setup machine info classes
    #---------------------------------------------------------------------------
    ats_args = install_ats_args()
    numNodes = options.numNodes
    timeLimit = options.timeLimit
    launch_cmd = ""
    blueOS = False
    # These are environment variables to suggest we are in an allocation already
    # NOTE: CI runs should already be in an allocation so the launch cmd is
    # unused in those cases
    inAllocVars = []

    if hostname:
        mac_args = []
        if any(x in hostname for x in toss_machine_names):
            numNodes = numNodes if numNodes else 2
            mac_args = [f"--numNodes {numNodes}"]
            timeLimit = timeLimit if timeLimit else 120
            inAllocVars = ["SLURM_JOB_NUM_NODES", "SLURM_NNODES"]
            launch_cmd = f"salloc --exclusive -N {numNodes} -t {timeLimit}"
        elif any(x in hostname for x in blueos_machine_names):
            blueOS = True
            numNodes = numNodes if numNodes else 2
            mac_args = ["--smpi_off", f"--numNodes {numNodes}"]
            inAllocVars = ["LSB_MAX_NUM_PROCESSORS"]
            timeLimit = timeLimit if timeLimit else 60
            launch_cmd = f"bsub -nnodes {numNodes} -Is -XF -core_isolation 2 -alloc_flags atsdisable -W {timeLimit}"
        if (options.ciRun):
            for i, j in ci_launch_flags.items():
                if (i in hostname):
                    launch_cmd += j + " "
        ats_args.extend(mac_args)

    #---------------------------------------------------------------------------
    # Launch ATS
    #---------------------------------------------------------------------------
    # If doing a CI run, set some more options
    if (options.ciRun):
        if ("--logs" not in unknown_options):
            ats_args.append(f"--logs {test_log_name}")
            log_name = test_log_name
        else:
            log_name_indx = unknown_options.index("--logs") + 1
            log_name = unknown_options[log_name_indx]
        ats_args.append("--continueFreq=15")
        # Pass flag to tell tests this is a CI run
        ats_args.append("--glue='cirun=True'")
    if (options.threads):
        ats_args.append(f"--glue='threads={options.threads}'")
    ats_args.append(f"""--glue='benchmark_dir="{benchmark_dir}"'""")
    ats_args.append("--glue='independent=True'")
    ats_args = " ".join(str(x) for x in ats_args)
    other_args = " ".join(str(x) for x in unknown_options)
    cmd = f"{ats_exe} -e {spheral_exe} {ats_args} {other_args}"
    # Check if are already in an allocation
    inAlloc = any(e in list(os.environ.keys()) for e in inAllocVars)
    # If already in allocation, do not do a launch
    if inAlloc:
        run_command = cmd
    else:
        if blueOS:
            # Launches using Bsub requires quoting the whole command
            # This causes issues for the glue='benchmark_dir... line
            # unless we escape the characters
            run_command = f'{launch_cmd} "{cmd}"'
            run_command = run_command.replace('="', '=\\"')
            run_command = run_command.replace('"\'', '\\"\'')
        else:
            run_command = f"{launch_cmd}{cmd}"
    print(f"\nRunning: {run_command}\n")
    if (options.ciRun):
        run_and_report(run_command, log_name, 0)
    else:
        try:
            subprocess.run(run_command, shell=True, check=True, text=True)
        except Exception as e:
            print(e)

if __name__ == "__main__":
    main()
