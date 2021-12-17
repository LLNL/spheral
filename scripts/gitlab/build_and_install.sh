#!/usr/bin/env bash

set -o errexit
set -o nounset

option=${1:-""}
hostname="$(hostname)"
project_dir="$(pwd)"

build_root=${BUILD_ROOT:-""}
spheral_spack_prefix=${SPHERAL_SPACK_PREFIX:-""}
mirror_dir=${SPHERAL_MIRROR_DIR:-""}
upstream_dir=${SPHERAL_UPSTREAM_DIR:-""}

hostconfig=${HOST_CONFIG:-""}
spec=${SPEC:-""}
job_unique_id=${CI_JOB_ID:-""}

sys_type=${SYS_TYPE:-""}
py_env_path=${PYTHON_ENVIRONMENT_PATH:-""}

lc_modules=${LC_MODULES:-""}

# Dependencies
date
if [[ "${option}" != "--build-only" && "${option}" != "--test-only" ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~~~~~ Building Dependencies"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    if [[ -z ${spec} ]]
    then
        echo "SPEC is undefined, aborting..."
        exit 1
    fi

    #prefix_opt=""
    #prefix="${project_dir}/../${hostname}"
    #if [[ -z ${job_unique_id} ]]; then
    #  job_unique_id=manual_job_$(date +%s)
    #  while [[ -d ${prefix}/${job_unique_id} ]] ; do
    #      sleep 1
    #      job_unique_id=manual_job_$(date +%s)
    #  done
    #fi
    #prefix="${prefix}/${job_unique_id}"

    echo "** Setting up spack"
    if [[ -z ${spheral_spack_prefix} ]]
    then
        spheral_spack_prefix="${project_dir}/../spheral-spack-tpls"
    fi

    mkdir -p ${spheral_spack_prefix}
    prefix_opt="--prefix=${spheral_spack_prefix}"


    echo "** Setting up upstream"
    upstream_opt=""
    if [[ ! -z ${upstream_dir} ]]
    then
        echo "** SPHERAL_UPSTREAM_DIR defined : Adding upstream ${upstream_dir}"
        upstream_opt="--upstream=${upstream_dir}"
    fi
    
    python3 scripts/uberenv/uberenv.py --setup-only ${prefix_opt} ${upstream_opt}



    echo "Activating spack instance..."
    . ${spheral_spack_prefix}/spack/share/spack/setup-env.sh


    echo "** Setting up mirror"
    if [[ ! -z ${mirror_dir} ]]
    then
        echo "** SPHERAL_MIRROR_DIR defined : Adding mirror ${mirror_dir}"
    else
        mirror_dir="/usr/WS2/davis291/SPHERAL/spheral-tpl/mirror"
        echo "** SPHERAL_MIRROR_DIR NOT defined : Adding mirror ${mirror_dir}"
    fi
    if spack mirror list | grep spheral-tpl; then
      echo "Spheral-tpl mirror found, removing..."
      spack mirror rm spheral-tpl
    fi
    echo "Adding Spheral-tpl mirror..."
    spack mirror add spheral-tpl ${mirror_dir}
    spack gpg trust `find ${mirror_dir} -name "*.pub"`


    echo "** Building TPL's and Host Config"
    echo "spack dev-build spheral@develop%${spec} hostconfig"
    #spack spec -I spheral@develop${spec}
    spack dev-build --quiet -d ${project_dir} -u initconfig spheral@develop%${spec} 2>&1 | tee -a dev-build-out.txt
fi
date


# Host config file
if [[ -z ${hostconfig} ]]
then
    # If no host config file was provided, we assume it was generated.
    # This means we are looking of a unique one in project dir.
    generated_hostconfig=${hostname//[0-9]/}-${sys_type}-${spec}
    generated_hostconfig_path=${project_dir}/${generated_hostconfig}.cmake
    echo "Searching for ${generated_hostconfig_path}..."
    
    hostconfigs=( $( ls "${project_dir}/"*.cmake ) )

    if [[ -f ${generated_hostconfig_path} ]]
    then
        hostconfig_path=${generated_hostconfig_path}
        hostconfig=${generated_hostconfig}
        echo "Found host config file: ${hostconfig_path}"
        echo "Host-config handle: ${hostconfig}"
        
        if [[ "${option}" != "--build_only" && ! -z "lc_modules" ]]
        then
            echo "#--------------------" >> ${hostconfig_path}
            echo "# LC_MODULES:${lc_modules}" >> ${hostconfig_path}
            echo "#--------------------" >> ${hostconfig_path}
        fi

    elif [[ ${#hostconfigs[@]} == 1 ]]
    then
        echo "Searching for unique host-config in ${project_dir}..."
        hostconfig_path=${hostconfigs[0]}
        echo "Found host config file: ${hostconfig_path}"
        hostconfig=$(echo ${hostconfig_path} | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
        echo "Host-config handle: ${hostconfig}"

    elif [[ ${#hostconfigs[@]} == 0 ]]
    then
        echo "No result for: ${project_dir}/*.cmake"
        echo "Spack generated host-config not found."
        exit 1

    else
        echo "More than one result for: ${project_dir}/*.cmake"
        echo "${hostconfigs[@]}"
        echo "Please specify one with HOST_CONFIG variable"
        exit 1
    fi
else
    # Using provided host-config file.
    hostconfig_path="${project_dir}/host-configs/${hostconfig}"
fi

# Build Directory
if [[ -z ${build_root} ]]
then
    build_root=$(pwd)
fi

build_dir="${build_root}/build_${hostconfig}/build"
install_dir="${build_root}/build_${hostconfig}/install"

cmake_exe=`grep 'CMake executable' ${hostconfig_path} | cut -d ':' -f 2 | xargs`


# Build
if [[ "${option}" != "--deps-only" && "${option}" != "--test-only" ]]
then
    date
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~ Host-config: ${hostconfig_path}"
    echo "~ Build Dir:   ${build_dir}"
    echo "~ Install Dir: ${install_dir}"
    echo "~ Project Dir: ${project_dir}"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo ""

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "~~~~~ Building Spheral"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    # If building, then delete everything first
    rm -rf ${build_dir} 2>/dev/null
    rm -rf ${install_dir} 2>/dev/null
    mkdir -p ${install_dir}
    mkdir -p ${build_dir} && cd ${build_dir}

    # Load modules for LC
    if [[ ! -z ${lc_modules} ]]
    then
        ml load ${lc_modules}
    else
        echo "Warning: No LC_MODULES set, ensure appropriate compilers are in path or you may experience incorrect builds!"
    fi

    date
    $cmake_exe \
      -C ${hostconfig_path} \
      -DCMAKE_INSTALL_PREFIX=${install_dir} \
      ${project_dir}
    if ! $cmake_exe --build . -j 48 --target install; then
      echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      echo "Compilation failed, running make VERBOSE=1"
      echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      $cmake_exe --build . --verbose --target install -j 1
    else
      ${install_dir}/spheral -c "import Spheral"
      #${install_dir}/spheral-atstest --atsExe /usr/gapps/ats/${sys_type}/7.0.5/bin/ats ${install_dir}/tests/integration.ats
      #/usr/gapps/ats/${sys_type}/7.0.5/bin/ats -e ${install_dir}/spheral ${install_dir}/tests/integration.ats
    fi
    date
fi
date
