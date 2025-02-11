#!/usr/bin/env python3

import argparse, os, sys, re, shutil

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from spheralutils import sexe

def get_config_dir(base_dir):
    "Return directory containing the repo.yaml file for a base dir"
    return os.path.join(base_dir, "scripts/spack")

# Spack instance info
os.environ["SPACK_DISABLE_LOCAL_CONFIG"] = "true"
default_spack_dir = os.path.join(os.getcwd(), "../spheral-spack-tpls")
default_spack_url = "https://github.com/spack/spack.git"
spack_commit = "5fe93fee1eec46a0750bd340198bffcb92ff9eec"
# Current repo (either LLNLSpheral or Spheral)
package_name = "spheral"

base_dir = os.getcwd()
package_dirs = {"spheral": base_dir}

# Find if this repo is LLNLSpheral by checking the submodule list
git_mod_cmd = "git config --file .gitmodules --name-only --get-regexp path$"
git_mod_out = sexe(git_mod_cmd, ret_output=True, echo=False)
if "spheral" in git_mod_out:
    package_name = "llnlspheral"
    package_dirs["spheral"] = os.path.join(base_dir, "spheral")
    package_dirs.update({"llnlspheral": base_dir})
print(f"Managing {package_name} TPLs")

dev_specs = {"gcc": "%gcc+mpi", "clang": "%clang+mpi", "gcc~mpi": "%gcc~mpi",
             "rocmcc": "+mpi~rocm", "rocmcc+rocm": "+mpi+rocm amdgpu_target=gfx942",
             "rocmcc+rocm~mpi": "~mpi+rocm amdgpu_target=gfx942"}

class SpheralTPL:
    def parse_args(self):
        parser = argparse.ArgumentParser()
        group = parser.add_mutually_exclusive_group()
        group.add_argument("--spec", type=str, default=None,
                           help="Install TPLs and create host config file for a given spec.\n"+\
                           "If no spec if given, TPLs for all env specs are installed and "+\
                           "no host config file is created.")
        group.add_argument("--dev-spec", type=str, default=None,
                           help="For use by developers to ease frequent TPL building",
                           choices=list(dev_specs.keys()))
        parser.add_argument("--show-spec", action="store_true",
                            help="Run spack spec -IL <spec>")
        parser.add_argument("--add-spec", action="store_true",
                            help="Set this flag to add the --spec to the environment.")
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
                            help="For use by the CI only. Must set a --spec.")
        parser.add_argument("--dry-run", action="store_true",
                            help="Use to do everything but actually install. For testing purposes.")
        parser.add_argument("--id", type=str, default=None,
                            help="ID string to postfix an initconfig file.")

        self.args = parser.parse_args()

        if (not self.args.spec and self.args.ci_run):
            raise Exception("Must specify a --spec if doing --ci-run")

        if (self.args.dev_spec):
            self.args.spec = package_name + dev_specs[self.args.dev_spec]

        if (self.args.spec and not self.args.spec.startswith(package_name)):
            raise Exception(f"--spec must start with {package_name}")
        if (self.args.spec):
            print(f"Installing {self.args.spec}")

    def add_spack_paths(self, spack_dir):
        "Append spack path to system to use spack python modules"
        spack_path = os.path.join(spack_dir, "lib", "spack")
        sys.path.append(spack_path)
        spack_external_path = os.path.join(spack_path, "external")
        sys.path.append(spack_external_path)
        sys.path.append(os.path.join(spack_external_path, "_vendoring"))
        global spack, SpackCommand
        try:
            import spack
            from spack.main import SpackCommand
            spack = spack
        except ImportError as e:
            raise ImportError("Failed to import Spack python module") from e

    def clone_spack(self):
        "Clone Spack and add paths to use spack python"
        tpl_root = self.args.spack_dir
        if (not os.path.exists(tpl_root)):
            os.mkdir(tpl_root)
        spack_dir = os.path.join(tpl_root, "spack")
        if (not self.args.skip_init):
            if (not os.path.exists(spack_dir)):
                sexe(f"git init {spack_dir}")
                sexe(f"git -C {spack_dir} remote add origin {self.args.spack_url}")
                sexe(f"git -C {spack_dir} fetch --depth=2 origin {spack_commit}")
                sexe(f"git -C {spack_dir} checkout FETCH_HEAD")
            else:
                # Check commit hash of Spack repo
                cur_hash = sexe(f"git -C {spack_dir} rev-parse HEAD", ret_output=True, echo=False).strip()
                if (cur_hash != spack_commit):
                    sexe(f"git -C {spack_dir} fetch --depth=2 origin {spack_commit}")
                    sexe(f"git -C {spack_dir} checkout FETCH_HEAD")
        self.add_spack_paths(spack_dir)
        if self.args.debug:
            sexe(f"git -C {spack_dir} clean -df")

    def custom_spack_env(self, env_dir, env_name):
        "Use/create a custom Spack environment"
        from spack import environment
        if (not self.args.spec):
            raise Exception("Must supply a --spec for a custom environment")

        cur_env_dir = os.path.join(env_dir, env_name)
        if (not os.path.exists(os.path.join(cur_env_dir, "spack.yaml"))):
            # Create a new environment
            env_cmd = SpackCommand("env")
            env_cmd("create", "--without-view", "-d", cur_env_dir)

        self.spack_env = environment.Environment(cur_env_dir)
        environment.activate(self.spack_env)
        # Get all the Spack commands
        repo_cmd = SpackCommand("repo")
        dev_cmd = SpackCommand("develop")
        comp_cmd = SpackCommand("compiler")
        ext_cmd = SpackCommand("external")

        # Add the repos and develop paths to the spack environment
        cur_repos = repo_cmd("list") # spack repo list
        for package, path in package_dirs.items():
            if (package+" " not in cur_repos):
                repo_path = os.path.abspath(get_config_dir(path))
                repo_cmd("add", f"{repo_path}") # spack repo add <repo_path>
            dev_path = os.path.abspath(path)
            dev_cmd("-p", dev_path, f"{package}@=develop") # spack develop <package>@=develop

        # Find external packages and compilers
        # List of packages to find externally
        ext_packages = ["git", "pkg-config", "autoconf", "automake"]
        if (not spack.spec.Spec(self.args.spec).satisfies("~mpi")):
            ext_packages.append("mpich")
        comp_cmd("find") # spack compiler find
        # Ignore any packages that are already found
        cur_packages = spack.config.get("packages")
        for i in ext_packages:
            if (i in cur_packages):
                ext_packages.remove(i)
        ext_cmd("find", *ext_packages)
        cur_packages = spack.config.get("packages")
        # TODO: Add logic to inform user when packages arent found
        # to encourage them to potentially add their own paths
        # Always add the spec for a custom environment
        self.args.add_spec = True

    def remove_upstream(self, env_dir):
        "Modify the spack.yaml to remove the upstream"
        # Copy original file
        from spack.util import spack_yaml
        env_file = os.path.join(env_dir, "spack.yaml")
        # TODO: Currently, Spack has no other way to
        # to remove an include: line from an environment
        # than to directly change the spack.yaml file

        # Load the spack.yaml file
        with open(env_file) as ff:
            try:
                loader = spack_yaml.load(ff)
            except SpackYAMLError as exception:
                print(exception)

        modded_file = False
        # Remove upstream.yaml or upstream entry
        if ("upstreams" in loader["spack"]):
            del loader["spack"]["upstreams"]
            modded_file = True
        if ("include" in loader["spack"]):
            for i, x in enumerate(loader["spack"]["include"]):
                if ("upstreams.yaml" in x):
                    del loader["spack"]["include"][i]
                    modded_file = True

        # Copy spack.yaml to origspack.yaml and overwrite spack.yaml
        # with upstreams removed
        if (modded_file):
            shutil.copyfile(env_file, os.path.join(env_dir, "origspack.yaml"))
            with open(env_file, 'w') as ff:
                spack_yaml.dump(loader, ff)

    def activate_spack_env(self):
        "Activates a Spack environment or creates and activates one when necessary"
        env_dir = os.path.join(get_config_dir(base_dir), "environments")
        # Check if we are on an LC machine and the environment exists
        default_env = os.getenv("SYS_TYPE")
        if default_env and os.path.exists(os.path.join(env_dir, default_env)):
            # For LC systems
            cur_env_dir = os.path.join(env_dir, default_env)
            print(f"Activating Spack environment in {cur_env_dir}")
            if self.args.no_upstream:
                self.remove_upstream(cur_env_dir)
            from spack import environment
            self.spack_env = environment.Environment(cur_env_dir)
            environment.activate(self.spack_env)
        else:
            # Otherwise, check if environment has been created
            arch_cmd = SpackCommand("arch")
            env_name = arch_cmd().strip()
            self.custom_spack_env(env_dir, env_name)

    def concretize_spec(self):
        "Concretize the spec"
        self.spack_spec = spack.spec.Spec(self.args.spec)
        if (self.args.add_spec):
            add_cmd = SpackCommand("add")
            add_cmd(self.args.spec)
        print("Concretizing environment")
        conc_cmd = SpackCommand("concretize")
        conc_cmd("-U")
        matches = self.spack_env.matching_spec(self.spack_spec)
        if (not matches):
            raise Exception(f"{self.args.spec} not found in current "+\
                            "environment. Rerun with --add-spec to add it.")
        self.spack_spec = matches
        print(f"Found matching root spec for {self.args.spec}")

    def install_and_config(self):
        "Install TPLs and create host config file for given spec"
        spec = self.args.spec
        # Load the spack package recipe python class
        if (package_name == "llnlspheral"):
            from spack.pkg.llnlspheral.llnlspheral import Llnlspheral
            spack_spheral = Llnlspheral(self.spack_spec)
        else:
            from spack.pkg.spheral.spheral import Spheral
            spack_spheral = Spheral(self.spack_spec)

        # Get host config file name from spack package recipe
        host_config_file = spack_spheral.cache_name
        # If using --id, preserve original host config file
        mod_host_config = False
        if (self.args.id and os.path.exists(host_config_file)):
            # Avoid overwriting existing host config file
            shutil.copyfile(host_config_file, "orig"+host_config_file)
            mod_host_config = True

        if (self.args.show_spec or self.args.ci_run):
            spec_cmd = SpackCommand("spec")
            print(f"Running spack spec -IL {spec}")
            spec_cmd("-IL", spec)
        if (not self.args.dry_run):
            install_cmd = SpackCommand("install")
            print(f"Running spack -u initconfig {spec}")
            install_cmd("-u", "initconfig", spec)
            print(f"Created {host_config_file}")

        if (self.args.ci_run):
            shutil.copyfile(host_config_file, "gitlab.cmake")

        if (mod_host_config):
            # Apply --id and bring back original host config file
            new_name = host_config_file.replace(".cmake", f"{self.args.id}.cmake")
            os.rename(host_config_file, new_name)
            os.rename("orig"+host_config_file, host_config_file)

    def __init__(self):
        self.parse_args()
        self.clone_spack()
        if (self.args.init_only):
            return
        self.activate_spack_env()
        if (self.args.spec):
            # If --spec is given, install TPLs and create host config file
            self.concretize_spec()
            self.install_and_config()
        else:
            # Concretize the current environment
            print("Concretizing environment")
            conc_cmd = SpackCommand("concretize")
            conc_cmd("-U")
            # No spec is given, install TPLs for all env specs
            install_cmd = SpackCommand("install")
            print(f"Running spack install --only dependencies {package_name}")
            install_cmd("--only", "dependencies", package_name)

        # Undo any file changes we made to spack.yaml
        orig_file = os.path.join(self.spack_env.path, "origspack.yaml")
        if (self.args.no_upstream and os.path.exists(orig_file)):
            # Revert env file if it was modified
            os.rename(orig_file, os.path.join(self.spack_env.path, "spack.yaml"))

if __name__=="__main__":
    spheral_tpl = SpheralTPL()
