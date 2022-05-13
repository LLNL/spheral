#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import json

#------------------------------------------------------------------------------

project_dir=os.path.abspath(os.path.join(os.path.realpath(__file__), "../../../"))

default_spheral_spack_dir=os.path.join(os.getcwd(), "../spheral-spack-tpls")
upstream_dir="/usr/WS2/wciuser/Spheral/spheral-spack-tpls/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_p/"

uberenv_path = os.path.join(project_dir, "scripts/devtools/uberenv/uberenv.py")
uberenv_project_json = os.path.join(os.getcwd(), ".uberenv_config.json")

#------------------------------------------------------------------------------

def parse_args():
  parser = argparse.ArgumentParser()

  # Spec args
  parser.add_argument('--spec', type=str, default="",
      help='Spack spec to use.')
  parser.add_argument('--spec-list', type=str, default="",
      help='JSON file with a list of specs to build for, this will override --spec.')

  # Spack setup
  parser.add_argument('--spheral-spack-dir', type=str, default=default_spheral_spack_dir,
      help='Dir of spack instance to handle tpls. This can point at \
            an existing spack dir or an empty on to create a new spack instance.')

  parser.add_argument('--no-clean', type=bool, default=False,
      help='Do not clean spack generated log files.')

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


# Parse the json formatted spec list...
def parse_spec_list(file_path):
  with open(file_path) as f:
    spec_list = json.loads(f.read())["specs"][os.environ.get("SYS_TYPE")]
  return spec_list


#------------------------------
# Dependencies
#------------------------------
def build_deps(args):
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("~~~~~ Building Dependencies")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("")
  print("{0}".format(project_dir))


  # Figure out what specs this script is building TPLs for.
  spec_list=[]
  if args.spec_list:
    spec_list = parse_spec_list(args.spec_list)
  elif args.spec:
    spec_list.append(args.spec)
  else:
    print("ERROR: Please define --spec or --spec-list, aborting...")
    sys.exit(1) 
  for s in spec_list:
    print("** SPEC : {0}".format(s))

  spack_config_dir_opt=""
  if "SYS_TYPE" not in os.environ.keys():
    spack_config_dir_opt="--spack-config-dir={0}".format(os.path.join(project_dir, "scripts/spack/configs/x86_64"))

  spack_upstream_opt=""
  if os.path.isdir(upstream_dir):
    spack_upstream_opt="--upstream {0}".format(upstream_dir)


  # We use uberenv to set up our spack instance with our respective package.yaml files
  # config.yaml and custom spack packages recipes.
  print("** Running uberenv...")
  prefix_opt="--prefix=" + args.spheral_spack_dir
  uberenv_project_json_opt="--project-json={0}".format(uberenv_project_json)

  print("** Spheral Spack Dir : {0}".format(args.spheral_spack_dir))
  sexe("python3 {0} --setup-only {1} {2} {3} {4}".format(uberenv_path, prefix_opt, uberenv_project_json_opt, spack_config_dir_opt, spack_upstream_opt), echo=True)

  # We just want to use the spac instance directly to generate our TPLs, we don't want
  # to have the spack instance take over our environment.
  os.environ["SPACK_DISABLE_LOCAL_CONFIG"] = "1"
  spack_cmd=os.path.join(args.spheral_spack_dir, "spack/bin/spack")

  with open(uberenv_project_json) as f:
    package_name=json.loads(f.read())["package_name"]

  # Loop through the specs we want TPLs for and build/install/get them as necessary.
  for s in spec_list:
    print("** Building TPL's and generating host-config for {0}%{1} ...".format(package_name,s))
    os.environ["SPEC"] = s
    if sexe("{0} spec -I {1}@develop%{2}".format(spack_cmd, package_name, s), echo=True) : sys.exit(1)
    if sexe("{0} dev-build --deprecated --quiet -d {1} -u initconfig {2}@develop%{3} 2>&1 | tee -a \"dev-build-{3}-out.txt\"".format(spack_cmd, os.getcwd(), package_name, s), echo=True) : sys.exit(1)

  if not args.no_clean:
    sexe("rm dev-build-* spack-build-* spack-configure-args.txt")

#------------------------------------------------------------------------------

def main():
  args = parse_args()
  build_deps(args)


if __name__ == "__main__":
  main()
