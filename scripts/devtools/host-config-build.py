#!/usr/bin/env python3

import os
import sys
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from spheralutils import sexe

source_dir=os.getcwd()

def parse_args():
  parser = argparse.ArgumentParser()

  parser.add_argument('--host-config', type=str, default="", required=True,
      help='host-config path.')

  parser.add_argument('-s', '--source-dir', type=str, default=source_dir,
      help='Location of spheral source directory.')

  parser.add_argument('-i', '--install-dir', type=str, default="",
      help='Location of spheral source directory.')

  parser.add_argument('--build-dir', type=str, default="",
      help='Name of build directory.')

  parser.add_argument('--no-clean', action='store_true',
      help='Do not delete build and install locations.')

  # Control stage flow
  parser.add_argument('--build', action='store_true',
      help='Run make -j install after configuring build dirs.')

  parser.add_argument('--lc-modules', type=str, default="",
      help='LC Modules to use during build, install and smoke test. This is not used if --build is not enabled.')

  parser.add_argument('-D', action='append', default=[])

  return parser.parse_args()

def main():
  args = parse_args()
  print(args)

  hostconfig=(args.host_config).split("/")[-1].split(".cmake")[0]
  if os.path.isabs(args.host_config):
    hostconfig_path=args.host_config
  else:
    hostconfig_path=os.path.abspath(args.host_config)


  # Set up our directory structure paths.
  if not args.build_dir:
    build_dir="{0}/build_{1}".format(source_dir, hostconfig)
  else:
    build_dir="{0}/{1}".format(source_dir, args.build_dir)
  if not args.install_dir:
    install_dir="{0}/install".format(build_dir)
  else:
    install_dir=args.install_dir
  build_dir=build_dir+"/build"
  # Pull the cmake command to use out of our host config.
  cmake_cmd=sexe("grep 'CMake executable' \"{0}\"".format(hostconfig_path), ret_output=True, echo=False).split()[-1]

  cmake_extra_args=""
  if args.D and args.D != ['']:
    cmake_extra_args="-D"+" -D".join(args.D)
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("~ Host-config:      {0}".format(hostconfig_path))
  print("~ Build Dir:        {0}".format(build_dir))
  print("~ Install Dir:      {0}".format(install_dir))
  print("~ Source Dir:       {0}".format(source_dir))
  print("~ Cmake cmd:        {0}".format(cmake_cmd))
  print("~ Extra CMake Args: {0}".format(cmake_extra_args))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("")

  # Clean our build and install dirs...
  if not args.no_clean:
      sexe("rm -rf \"{0}\" 2>/dev/null".format(build_dir),echo=True)
      sexe("rm -rf \"{0}\" 2>/dev/null".format(install_dir),echo=True)
  sexe("mkdir -p \"{0}\"".format(build_dir),echo=True)
  sexe("mkdir -p \"{0}\"".format(install_dir),echo=True)

  # Move to the build directory.
  os.chdir(build_dir)
  
  # We need a module command to add to our environment on each python
  # subprocess call otherwise we might get a bad build depending on 
  # the default compiler setup for the system.
  ml_cmd=""
  ml_cmd_v=""
  if not args.lc_modules:
    print("Warning: No LC_MODULES set, ensure appropriate compilers are in path or you may experience incorrect builds!")
  else:
    ml_cmd_v = "module load {0} &&".format(args.lc_modules)
    ml_cmd = "module load {0} 2>/dev/null &&".format(args.lc_modules)

  # Run our CMake config step.
  sexe("{0} {1} -C {2} -DCMAKE_INSTALL_PREFIX={3} {4} {5}".format(ml_cmd_v, cmake_cmd, hostconfig_path, install_dir, cmake_extra_args, source_dir))

  # Build and install Spheral

  if args.build:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~ Building Spheral")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    sexe("{0} {1} --build . -j 48 --target install".format(ml_cmd, cmake_cmd), echo=True, ret_output=False)

if __name__ == "__main__":
  main()
