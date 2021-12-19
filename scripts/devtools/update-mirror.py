#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import json


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--spheral-spack-dir', type=str, required=True,
      help='Dir of spack instance where TPLs are installed.')

  parser.add_argument('--mirror-dir', type=str, required=True,
      help='Dir of mirror to be used when --use-mirror is enabled.')

  #parser.add_argumtnt('--generate-new-key', type=bool, default=False,
  #    help='Should spack generate a new key, only use this if creating a new mirror.')

  parser.add_argument('--secret-key-dir', type=str, required=True,
      help='Dir where the secret keys are stored.')


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

def update_mirror(args):
  os.environ["SPACK_DISABLE_LOCAL_CONFIG"] = "1"
  spack_cmd=os.path.join(args.spheral_spack_dir, "spack/bin/spack")

  if sexe("{0} mirror list | grep spheral-tpl".format(spack_cmd), echo=True, ret_output=True):
    print("** Spheral-tpl mirror found, removing...")
    sexe("{0} mirror rm spheral-tpl".format(spack_cmd), echo=True)


  spack_build_cache_dir=os.path.join(args.mirror_dir, "build_cache")
  spack_gpg_dir=os.path.join(args.spheral_spack_dir, "spack/opt/spack/gpg")
  spack_gpg_backup_dir=os.path.join(args.spheral_spack_dir, "spack/opt/spack/gpg-backup")

  sexe("mkdir -p \"{0}\"".format(args.mirror_dir), echo=True)
  sexe("mkdir -p \"{0}\"".format(spack_build_cache_dir), echo=True)


  sexe("mkdir -p \"{0}\"".format(spack_gpg_backup_dir), echo=True)

  sexe("mv {0}/* {1}".format(spack_gpg_dir, spack_gpg_backup_dir), echo=True)
  sexe("cp {0}* {1}".format(args.secret_key_dir, spack_gpg_dir), echo=True)


  # Create Source Mirror
  #sexe("{0} mirror create --directory {1}".format(spack_cmd, args.mirror_dir), echo=True)
  #sexe("chmod -R g+rws {0}".format(args.mirror_dir), echo=True)

  sexe("{0} mirror add spheral-tpl {1}".format(spack_cmd, args.mirror_dir), echo=True)
  sexe("{0} buildcache keys --install --trust".format(spack_cmd), echo=True)

  sexe("for ii in $({0} find --format \"yyy {{name}} {{version}} /{{hash}}\" | grep -v -E \"^(develop^master)\" | grep -v spheral | grep \"yyy\" | cut -f4 -d \" \" ); do {0} buildcache create --allow-root --force -d {1} --only=package $ii; done".format(spack_cmd, args.mirror_dir), echo=True)
  sexe("chmod -R g+rws {0}".format(spack_build_cache_dir), echo=True)

  

def main():
  args = parse_args()
  update_mirror(args)


if __name__ == "__main__":
  main()
