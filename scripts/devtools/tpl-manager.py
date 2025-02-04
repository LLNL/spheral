#!/usr/bin/env python3

import argparse, os, sys, re, yaml, shutil

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from spheralutils import sexe

# Spack instance info
os.environ["SPACK_DISABLE_LOCAL_CONFIG"] = "true"
default_spack_dir = os.path.join(os.getcwd(), "../spheral-spack-tpls")
default_spack_url = "https://github.com/spack/spack.git"
spack_commit = "5fe93fee1eec46a0750bd340198bffcb92ff9eec"
# Current repo (either LLNLSpheral or Spheral)
base_dir = os.getcwd()
spheral_config_dir = os.path.join(base_dir, "scripts/spack")
package_name = "spheral"
# Find if this repo is LLNLSpheral by checking the submodule list
git_mod_cmd = "git config --file .gitmodules --name-only --get-regexp path$"
git_mod_out = sexe(git_mod_cmd, ret_output=True)
if "spheral" in git_mod_out:
    package_name = "llnlspheral"

# Default name of spack environment on LC systems
default_spack_env = os.getenv("SYS_TYPE")

class SpheralTPL:
    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--spec", type=str,
                            default=None, help="Spack spec to install.")
        parser.add_argument("--spack-dir", type=str,
                            default=default_spack_dir,
                            help="Directory to install Spack instance and a build directory.")
        parser.add_argument("--spack-url", type=str, default=default_spack_url,
                            help="URL to download spack.")
        parser.add_argument("--debug", action="store_true",
                            help="Set this flag if seeing odd Spack issues with tpl-manager.")
        parser.add_argument("--no-upstream", action="store_true",
                            help="Ignore upstream by temporarily modifying environment.")
        parser.add_argument("--init-only", action="store_true",
                            help="Download Spack but do not concretize or install.")
        parser.add_argument("--skip-init", action="store_true",
                            help="Skip downloading Spack repo.")
        parser.add_argument("--ci-run", action="store_true",
                            help="For use by the CI only.")
        parser.add_argument("--id", type=str, default=None,
                            help="ID string to postfix an initconfig file.")
        self.args = parser.parse_args()

    def clone_spack(self):
        "Clone Spack"
        tpl_root = self.args.spack_dir
        if (not os.path.exists(tpl_root)):
            os.mkdir(tpl_root)
        spack_dir = os.path.join(tpl_root, "spack")
        if (not self.args.skip_init):
            if (not os.path.exists(spack_dir)):
                sexe(f"git init {spack_dir}")
                sexe(f"git -C {spack_dir} remote add origin {self.args.spack_url}")
                sexe(f"git -C {spack_dir} config feature.manyFiles true")
            # Check commit hash of Spack repo
            cur_hash = sexe(f"git -C {spack_dir} rev-parse FETCH_HEAD", ret_output=True, echo=False).strip()
            if (cur_hash != spack_commit):
                sexe(f"git -C {spack_dir} fetch --depth=2 origin {spack_commit}")
                sexe(f"git -C {spack_dir} checkout FETCH_HEAD")
        self.spack_cmd = os.path.join(spack_dir, "bin/spack")
        if self.args.debug:
            sexe(f"git -C {spack_dir} clean -df")
            self.spack_cmd += " -d"

    def create_spack_env(self, env_dir, arch):
        "Create a new Spack environment from scratch for Spheral"
        os.chdir(env_dir)
        sexe(f"{self.spack_cmd} env create {arch} -d")
        self.spack_cmd += f" -e {os.path.join(env_dir, arch)}"
        find_cmd = "external find --non-buildable"
        spack_env_cmds = ["compiler find", find_cmd]
        mpi_packages = ["mvapich", "mvapich2", "spectrum-mpi", "cray-mpich"]
        spack_env_cmds.extend([f"{find_cmd} {x}" for x in mpi_packages])
        # TODO: Put additional external find commands for compilers and packages
        for i in spack_env_cmds:
            sexe(f"{self.spack_cmd} {i}")
        os.chdir(base_dir)

    def remove_upstream(self, cur_env):
        with open(self.env_file) as ff:
            try:
                loader = yaml.safe_load(ff)
            except yamlYAMLError as exception:
                print(exception)
        if ("upstreams" in loader["spack"]):
            # Copy original file to revert to
            self.env_file = os.path.join(cur_env, "spack.yaml")
            self.orig_env_file = os.path.join(cur_env, "origspack.yaml")
            shutil.copyfile(self.env_file, self.orig_env_file)
            del loader["spack"]["upstreams"]
            with open(spackyaml, 'w') as ff:
                yaml.dump(loader, ff)
        else:
            self.orig_env_file = None

    def activate_spack_env(self):
        "Activates a Spack environment by putting -e env_dir after all spack commands"
        # If we are on LC systems
        env_dir = os.path.join(spheral_config_dir, "environments")
        if default_spack_env:
            cur_env = os.path.join(env_dir, default_spack_env)
            if (not os.path.exists(cur_env)):
                raise Exception(f"{cur_env} does not exists")
            self.spack_cmd += f" -e {cur_env}"
            self.add_spec = False
        else:
            # Otherwise, check if environment has been created
            arch = sexe(f"{self.spack_cmd} arch", ret_output=True).strip()
            cur_env = os.path.join(env_dir, arch)
            if not os.path.exists(cur_env):
                # Create the environment
                self.create_spack_env(env_dir, arch)
            else:
                self.spack_cmd += f" -e {cur_env}"
            self.add_spec = True
        # To turn upstream off, copy the spack.yaml and modify it
        if self.args.no_upstream:
            self.remove_upstream(env_dir)

    def concretize(self):
        if self.args.debug:
            sexe(f"{self.spack_cmd} concretize -f")
        else:
            sexe(f"{self.spack_cmd} concretize")

    def install_tpls(self):
        spec = self.args.spec
        install_cmd = f"{self.spack_cmd} install -u initconfig")
        if (self.add_spec):
            install_cmd += " --add"
        if self.args.spec:
            install_cmd += f" {self.args.spec}"
        else:
            install_cmd += " spheral"
        sexe(install_cmd)

    def __init__(self):
        self.parse_args()
        if self.args.debug:
            os.rmdir("~/.spack")
        self.clone_spack()
        self.activate_spack_env()
        if (self.args.spec):
            if (not self.args.spec.startswith(package_name+"%")):
                self.args.spec = package_name + "%" + self.args.spec
            sexe(f"{self.spack_cmd} spec {self.args.spec}")
        # Name of initconfig file created from Spack
        host_config = f"{sys_type}-{self.args.spec.replace(' ', '_')}"
        host_config_file = host_config + ".cmake"
        if (not self.args.init_only):
            self.orig_host_config = None
            if (self.args.id and os.path.exists(host_config_file)):
                self.orig_host_config = "orig"+host_config_file
                shutil.copyfile(host_config_file, self.orig_host_config)
            self.concretize()
            self.install_tpls()
        # Undo any temporary file changes we made
        # for example to spack.yaml or host config file
        if (self.args.no_upstream and self.orig_env_file):
            # Revert env file if it was modified
            os.rename(self.orig_env_file, self.env_file)
        if (self.args.ci_run):
            shutil.copyfile(host_config, "gitlab.cmake")
        elif (self.orig_host_config):
            os.rename(host_config_file, host_config+self.args.id+".cmake")
            os.rename(self.orig_host_config, host_config_file)

if __name__=="__main__":
    spheral_tpl = SpheralTPL()
