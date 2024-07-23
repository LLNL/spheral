#!/usr/bin/env python3

from string import digits
import os
import sys
import argparse
import subprocess

#------------------------------------------------------------------------------

project_dir=os.path.abspath(os.path.join(os.path.realpath(__file__), "../../../"))
default_spheral_spack_dir=os.path.join(project_dir, "../spheral-spack-tpls")

tpl_manager_cmd=os.path.join(project_dir, "scripts/devtools/tpl-manager.py")
host_congfig_build_cmd=os.path.join(project_dir, "scripts/devtools/host-config-build.py")

host=os.environ.get("LCSCHEDCLUSTER")
sys_type=os.environ.get("SYS_TYPE")

#------------------------------------------------------------------------------

def parse_args():
  parser = argparse.ArgumentParser()

  # Control stage flow
  parser.add_argument('--tpls-only', action='store_true',
      help='Only build the tpls and generate the host-config.')
  parser.add_argument('--build-only', action='store_true',
      help='Only build the project from a spack generated host-config.')
  parser.add_argument('--host-config', type=str, default="",
      help='host-config path.')

  # Spec args
  parser.add_argument('--spec', type=str, default="", required=True,
      help='Spack spec to use.')
  parser.add_argument('--lc-modules', type=str, default="",
      help='LC Modules to use during build, install and smoke test.')

  parser.add_argument('--extra-cmake-args', type=str, default="",
      help='CMake args to pass to the build stage in the host-config-build script.')

  return parser.parse_args()


# Helper function for executing commands stolen from uberenv
def sexe(cmd,ret_output=False,echo=True):
    """ Helper for executing shell commands. """
    if echo:
      print("[exe: {0}]".format(cmd))
    p = subprocess.run(cmd, shell=True,
                       capture_output=ret_output,
                       check=True, text=True)
    if ret_output:
      if echo:
        print(p.stdout)
      return p.stdout

#------------------------------------------------------------------------------

def main():
  args = parse_args()

  # Build out our TPLs for this spec if we need to.
  if not args.build_only:
    print("** Building spec TPLs : {0}".format(args.spec))
    sexe("{0} --spec=\"{1}\"".format(tpl_manager_cmd, args.spec))
    #sexe("{0} --spheral-spack-dir spheral-spack-tpls --spec=\"{1}\"".format(tpl_manager_cmd, args.spec))

  # Get the host-config name and path.
  if not args.build_only and not args.host_config:
    hostconfig="{1}-{2}.cmake".format(host, sys_type, (args.spec).replace(" ","_"))
    sexe("cp {0} gitlab.cmake".format(hostconfig))
    hostconfig_path=os.path.join(os.getcwd(), "gitlab.cmake")
  else:
    hostconfig=(args.host_config).split("/")[-1].split(".cmake")[0]
    hostconfig_path=args.host_config
  print(hostconfig)

  if not args.tpls_only:
      if sexe("{0} --host-config=\"{1}\" --lc-modules=\"{2}\" --build {3}".format(host_congfig_build_cmd, hostconfig_path, args.lc_modules, args.extra_cmake_args)) : sys.exit(1)

if __name__ == "__main__":
  main()
