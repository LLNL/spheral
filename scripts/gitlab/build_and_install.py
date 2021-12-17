#!/usr/bin/env python

from string import digits
import os
import sys
import argparse
import subprocess

#------------------------------------------------------------------------------

project_dir=os.getcwd()
default_spheral_spack_dir=os.path.join(project_dir, "../spheral-spack-tpls")

tpl_manager_cmd=os.path.join(project_dir, "scripts/devtools/tpl-manager.py")

host=os.environ.get("HOSTNAME").translate(None, digits)
sys_type=os.environ.get("SYS_TYPE")

#------------------------------------------------------------------------------

def parse_args():
  parser = argparse.ArgumentParser()

  # Control stage flow
  parser.add_argument('--build-only', action='store_true',
      help='Only build the project from a spack generated host-config.')
  parser.add_argument('--host-config', type=str, default="",
      help='host-config path.')

  # Spec args
  parser.add_argument('--spec', type=str, default="", required=True,
      help='Spack spec to use.')
  parser.add_argument('--lc-modules', type=str, default="",
      help='LC Modules to use during build and install and smoke test.')

  return parser.parse_args()


# Helper function for executing commands stolen from uberenv
def sexe(cmd,ret_output=False,echo=False):
    """ Helper for executing shell commands. """
    if echo:
        print("[exe: {0}]".format(cmd))
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        out = p.communicate()[0]
        out = out.decode('utf8')
        return p.returncode,out
    else:
        return subprocess.call(cmd,shell=True)


#------------------------------------------------------------------------------

def main():
  args = parse_args()

  # Build out our TPLs for this spec if we need to.
  if not args.build_only:
    print("** Building spec TPLs : {0}".format(args.spec))
    sexe("{0} --spec={1}".format(tpl_manager_cmd, args.spec))

  # Get the host-config name and path.
  if not args.build_only and not args.host_config:
    hostconfig="{0}-{1}-{2}".format(host, sys_type, args.spec)
    hostconfig_path=os.path.join(project_dir, "{0}.cmake".format(hostconfig))
  else:
    hostconfig=(args.host_config).split("/")[-1].split(".cmake")[0]
    hostconfig_path=args.host_config
  print(hostconfig)

  # Set up our directory structure paths.
  build_dir="{0}/build_{1}/build".format(project_dir, hostconfig)
  install_dir="{0}/build_{1}/install".format(project_dir, hostconfig)

  # Pull the cmake command to use out of our host config.
  cmake_cmd=sexe("grep 'CMake executable' {0}".format(hostconfig_path), ret_output=True, echo=True)[1].split()[-1]
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("~ Host-config: {0}".format(hostconfig_path))
  print("~ Build Dir:   {0}".format(build_dir))
  print("~ Install Dir: {0}".format(install_dir))
  print("~ Project Dir: {0}".format(project_dir))
  print("~ Cmake cmd:   {0}".format(cmake_cmd))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("")

  # Clean our build and install dirs...
  sexe("rm -rf {0} 2>/dev/null".format(build_dir),echo=True)
  sexe("rm -rf {0} 2>/dev/null".format(install_dir),echo=True)
  sexe("mkdir -p {0}".format(build_dir),echo=True)
  sexe("mkdir -p {0}".format(install_dir),echo=True)

  # Move to the build directory.
  os.chdir(build_dir)
  
  # We need a module command to add to our environment on each python
  # subprocess call otherwise we might get a bad build depending on 
  # the default compiler setup for the system.
  ml_cmd=""
  if not args.lc_modules:
    print("Warning: No LC_MODULES set, ensure appropriate compilers are in path or you may experience incorrect builds!")
  else:
    ml_cmd = "module load {0} 2>/dev/null &&".format(args.lc_modules)


  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("~~~~~ Building Spheral")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  # Run our CMake config step.
  sexe("{0} {1} -C {2} -DCMAKE_INSTALL_PREFIX={3} {4}".format(ml_cmd, cmake_cmd, hostconfig_path, install_dir, project_dir))

  # Build and install Spheral
  build_result = sexe("{0} {1} --build . -j 48 --target install".format(ml_cmd, cmake_cmd), echo=True)

  # If our build or install failed, run again to get our first error.
  if build_result != 0:
    print(build_result)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Compilation failed, running make VERBOSE=1")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    sexe("{0} {1} --build . --verbose --target install -j 1".format(ml_cmd, cmake_cmd),ret_output=True, echo=True)
    sys.exit(1)

  # Try to import Spheral for a basic sanity test.
  smoke_test = sexe("{0} {1}/spheral -c \"import Spheral\"".format(ml_cmd, install_dir))
  if smoke_test != 0:
    sys.exit(1)



if __name__ == "__main__":
  main()
