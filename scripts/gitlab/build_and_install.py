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
host_congfig_build_cmd=os.path.join(project_dir, "scripts/devtools/host-config-build.py")

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
      help='LC Modules to use during build, install and smoke test.')

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
    sexe("{0} --spec=\"{1}\" --use-mirror True".format(tpl_manager_cmd, args.spec))

  # Get the host-config name and path.
  if not args.build_only and not args.host_config:
    hostconfig="{0}-{1}-{2}".format(host, sys_type, args.spec)
    hostconfig_path=os.path.join(project_dir, "{0}.cmake".format(hostconfig))
  else:
    hostconfig=(args.host_config).split("/")[-1].split(".cmake")[0]
    hostconfig_path=args.host_config
  print(hostconfig)

  sexe("{0} --host-config=\"{1}\" --lc-modules=\"{2}\" --build".format(host_congfig_build_cmd, hostconfig_path, args.lc_modules))


if __name__ == "__main__":
  main()