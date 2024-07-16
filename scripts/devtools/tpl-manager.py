#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import json

#------------------------------------------------------------------------------

project_dir=os.path.abspath(os.path.join(os.path.realpath(__file__), "../../../"))

default_spheral_spack_dir=os.path.join(os.getcwd(), "../spheral-spack-tpls")
default_upstream_dir="/usr/WS2/sduser/Spheral/spack_upstream/0.22/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_p"

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

  parser.add_argument('--upstream-dir', type=str, default=default_upstream_dir,
      help='Dir of upstream spack installation.')

  parser.add_argument('--no-upstream', action="store_true",
      help='Do not use an upstream spack instance.')

  parser.add_argument('--init-only', action="store_true",
      help='Only initialize the spack instance.')

  parser.add_argument('--spack-url', type=str, default="",
      help='URL of spack to use.')

  parser.add_argument('--spack-jobs', type=str,
      help='Optionally launch spack builds in parallel with given number of parallel jobs.')

  parser.add_argument('--no-clean', action='store_true',
      help='Do not clean spack generated log files.')

  parser.add_argument('--debug', action='store_true',
      help='Run spack commands with -d debug output.')

  parser.add_argument('--verbose', action='store_true',
      help='Run spack install with -v for tpl build output.')

  parser.add_argument('--no-spec', action='store_true',
      help='Skip output of the dependency graph.')

  parser.add_argument('--skip-init', action='store_true',
      help='Skip setting up and configuring Spack.')

  return parser.parse_args()


# Helper function for executing commands stolen from uberenv
def sexe(cmd,ret_output=False,echo=True):
    """ Helper for executing shell commands. """
    if echo:
      print("[exe: {0}]".format(cmd))

    # If we want to return the output as string a print to stdout
    # in real-time we need to let subprocess print as normal to 
    # PIPE and STDOUT. We then need to read it back ourselves and 
    # append to an ouput string of our own making. There is no way
    # to do this with subprocess currently.
    if ret_output:
      p = subprocess.Popen(cmd,
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           encoding='utf8')
      out = "";
      while True:
        realtime_output = p.stdout.readline()

        if realtime_output == '' and p.poll() is not None:
          break

        if realtime_output:
          print(realtime_output.strip(), flush=True)
          out += realtime_output

      if echo:
        print(out)
      return out
    
    # If we do not need to return the output as a string, run() 
    # will suffice.
    else:
      p = subprocess.run(cmd, shell=True,
                         check=True, text=True)


# Parse the json formatted spec list...
def parse_spec_list(file_path):
  with open(file_path) as f:
    spec_list = json.loads(f.read())["specs"][os.environ.get("SYS_TYPE")]
  return spec_list


#------------------------------
# Dependencies
#------------------------------
def build_spack(args):
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("~~~~~ Configuring Spack")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("")
  print("{0}".format(project_dir))

  spack_upstream_opt=""
  if os.path.isdir(args.upstream_dir) and not args.no_upstream:
    spack_upstream_opt="--upstream {0}".format(args.upstream_dir)

  uberenv_spack_url_opt=""
  if args.spack_url:
      uberenv_spack_url_opt="--spack-url {0}".format(args.spack_url)

  # We use uberenv to set up our spack instance with our respective package.yaml files
  # config.yaml and custom spack packages recipes.
  print("** Running uberenv...")

  prefix_opt="--prefix=" + args.spheral_spack_dir
  uberenv_project_json_opt="--project-json={0}".format(uberenv_project_json)

  print("** Spheral Spack Dir : {0}".format(args.spheral_spack_dir))

  # We just want to use the spac instance directly to generate our TPLs, we don't want
  # to have the spack instance take over our environment.
  os.environ["SPACK_DISABLE_LOCAL_CONFIG"] = "1"
  spack_cmd=os.path.join(args.spheral_spack_dir, "spack/bin/spack")

  spheral_config_dir="scripts/spack/configs/"
  spack_config_dir_opt=""
  if "SYS_TYPE" not in os.environ.keys():
    # We need to install spack without any configuration files so we can use
    # spack arch to determine the OS of the system and later to use spack find
    # for generating external package files on external systems.
    sexe("python3 {0} --setup-only {1} {2} {3} {4}".format(uberenv_path, prefix_opt, uberenv_project_json_opt, spack_upstream_opt, uberenv_spack_url_opt))

    spack_arch_os = sexe("{0} arch -o".format(spack_cmd), ret_output=True, echo=False).strip()
    print("INFO : Detected Operating System :{0}".format(spack_arch_os))

    spheral_config_dir += spack_arch_os

    spack_config_dir_opt="--spack-config-dir={0}".format(os.path.join(project_dir, spheral_config_dir))
  else:
    spheral_config_dir += os.environ["SYS_TYPE"]


  # Setup spack w/ Uberenv and the appropriate external package/compiler configs.
  sexe("python3 {0} --setup-only {1} {2} {3} {4} {5}".format(uberenv_path, prefix_opt, uberenv_project_json_opt, spack_config_dir_opt, spack_upstream_opt, uberenv_spack_url_opt))

  # Uberenv doesn't copy the concretizer.yaml options...
  if os.path.exists(spheral_config_dir+"/concretizer.yaml"):
    sexe("cp {0}/concretizer.yaml {1}".format(spheral_config_dir, os.path.join(args.spheral_spack_dir, "spack/etc/spack/defaults")))

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
  spack_cmd=os.path.join(args.spheral_spack_dir, "spack/bin/spack")
  
  # Optionally add a parallel job number for spack builds
  if args.spack_jobs:
    spack_cmd += " --jobs={0}".format(args.spack_jobs)

  # Add -d to spack command when requesting debug output
  if args.debug:
    spack_cmd += " -d"

  with open(uberenv_project_json) as f:
    package_name=json.loads(f.read())["package_name"]

  # Loop through the specs we want TPLs for and build/install/get them as necessary.
  if not args.init_only:
      for s in spec_list:
        print("** Building TPL's and generating host-config for {0}%{1} ...".format(package_name,s))
        os.environ["SPEC"] = s
        os.environ["LC_ALL"] = "en_US.UTF-8"

        if not args.no_spec:
          sexe("{0} spec --fresh -IL {1}@develop%{2} 2>&1 | tee -a \"spec-info-{2}-out.txt\"".format(spack_cmd, package_name, s))

        # Install only the dependencies for Spheral and create CMake configure file
        if "Error: " in sexe("{0} dev-build -q --fresh -u initconfig {1}@develop%{2} 2>&1 | tee -a \"tpl-build-{2}-out.txt\"".format(spack_cmd, package_name, s), ret_output=True): sys.exit(1)

      if not args.no_clean:
        sexe("rm -f spec-info-* tpl-build-* spack-build-* spack-configure-args.txt")

#------------------------------------------------------------------------------

def main():
  args = parse_args()
  if (not args.skip_init):
    build_spack(args)
  build_deps(args)


if __name__ == "__main__":
  main()
