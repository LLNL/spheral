#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import json

#------------------------------------------------------------------------------

project_dir=os.getcwd()
default_mirror_dir="/usr/WS2/davis291/gapps/Spheral/spheral-spack-tpls/mirror"
default_spheral_spack_dir=os.path.join(project_dir, "../spheral-spack-tpls")

#------------------------------------------------------------------------------

def parse_args():
  parser = argparse.ArgumentParser()

  # Spec args
  parser.add_argument('--spec', type=str, default="",
      help='Spack spec to use.')
  parser.add_argument('--spec-list', type=str, default="",
      help='JSON file with a list of specs to build for, this will override --spec.')

  # Mirrors
  parser.add_argument('--no-mirror', action='store_true',
      help='Use a mirror with the spack instancedt.')
  parser.add_argument('--mirror-dir', type=str, default=default_mirror_dir,
      help='Dir of mirror to be used when --use-mirror is enabled.')

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

  # We use uberenv to set up our spack instance with our respective package.yaml files
  # config.yaml and custom spack packages recipes.
  print("** Running uberenv...")
  prefix_opt="--prefix=" + args.spheral_spack_dir
  print("** Spheral Spack Dir : {0}".format(args.spheral_spack_dir))
  sexe("python3 scripts/devtools/uberenv/uberenv.py --setup-only {0} {1}".format(prefix_opt, spack_config_dir_opt))

  # We just want to use the spac instance directly to generate our TPLs, we don't want
  # to have the spack instance take over our environment.
  os.environ["SPACK_DISABLE_LOCAL_CONFIG"] = "1"
  spack_cmd=os.path.join(args.spheral_spack_dir, "spack/bin/spack")
  
  # Let's set up a mirror for our TPL builds, this will help in downloding tars for packages 
  # offline and for pulling in binaries of precompiled TPL specs.
  if not args.no_mirror:

    print("** Setting up mirror")
    if args.mirror_dir:
      print("** --mirro-dir defined : Adding mirror : {0}".format(args.mirror_dir))
    else:
      print("** --mirror-dir NOT defined : Adding default mirror : {0}".format(args.mirror_dir))

    if sexe("{0} mirror list | grep spheral-tpl".format(spack_cmd), ret_output=True):
      print("** Spheral-tpl mirror found, removing...")
      sexe("{0} mirror rm spheral-tpl".format(spack_cmd))

    print("** Adding Spheral-tpl mirror...")
    sexe("{0} mirror add spheral-tpl {1}".format(spack_cmd, args.mirror_dir))
    sexe("{0} gpg trust `find {1} -name \"*.pub\"`".format(spack_cmd, args.mirror_dir))


  # Loop through the specs we want TPLs for and build/install/get them as necessary.
  for s in spec_list:
    print("** Building TPL's and generating host-config for {0} ...".format(s))
    os.environ["SPEC"] = s
    sexe("{0} spec -I spheral@develop%{1}".format(spack_cmd, s))
    sexe("{0} dev-build --quiet -d {1} -u initconfig spheral@develop%{2} 2>&1 | tee -a \"dev-build-{2}-out.txt\"".format(spack_cmd, project_dir, s))

  if not args.no_clean:
    sexe("rm dev-build-* spack-build-* spack-configure-args.txt")

#------------------------------------------------------------------------------

def main():
  args = parse_args()
  build_deps(args)


if __name__ == "__main__":
  main()
