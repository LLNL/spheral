#!/usr/bin/env bash

set -o errexit
set -o nounset

option=${1:-""}
hostname="$(hostname)"
spec=${SPEC:-""}
job_unique_id=${CI_JOB_ID:-""}

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

    upstream_opt=""
    prefix_opt=""


    prefix="../${hostname}"
    if [[ -z ${job_unique_id} ]]; then
      job_unique_id=manual_job_$(date +%s)
      while [[ -d ${prefix}/${job_unique_id} ]] ; do
          sleep 1
          job_unique_id=manual_job_$(date +%s)
      done
    fi

    prefix="${prefix}/${job_unique_id}"
    spack_env_dir="${prefix}/spack-env"
    mirror_dir="/usr/WS2/davis291/SPHERAL/spheral-tpl/mirror"
    src_dir=$(pwd)

    #prefix="/usr/WS2/davis291/SPHERAL/spheral-pr/shm/rzgenie47/manual_job_1637216151"

    mkdir -p ${prefix}
    mkdir -p ${spack_env_dir}
    prefix_opt="--prefix=${prefix}"

    #upstream_opt="--upstream=/usr/WS2/davis291/SPHERAL/lc_uberenv_tpl2/spack/opt/spack"
    #upstream_opt="--upstream=/usr/workspace/wsb/davis291/SPHERAL/spheral-spack-tpls/uberenv_libs/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__"
    #fi

    #python3 scripts/uberenv/uberenv.py --spec="${spec}" --reuse=True ${upstream_opt} ${prefix_opt}
    python3 scripts/uberenv/uberenv.py --setup-only ${prefix_opt}

    echo "Activating spack instance..."
    . ${prefix}/spack/share/spack/setup-env.sh

    echo "Create and activate Spack environment..."
    cd ${spack_env_dir}
    spack env create -d .
    spack env activate .

    echo "Adding Spheral-tpl mirror..."
    spack mirror add spheral-tpl ${mirror_dir}
    spack gpg trust `find ${mirror_dir} -name "*.pub"`

    echo "spack dev-build spheral@develop${spec}"
    spack dev-build -d ${src_dir} spheral@develop${spec}

    echo "Activating spack instance..."
    . ${prefix}/spack/share/spack/setup-env.sh

    spheral_spack_hash=`spack find -L spheral | grep spheral | cut -d ' ' -f1`
    spheral_install_prefix=`spack find -p spheral | grep spheral | cut -d ' ' -f3`

    echo "Loading spheral ${spheral_spack_hash} ..."
    spack load /${spheral_spack_hash}

    cd ${spheral_install_prefix}

    echo "Spheral import test..."
    ./spheral -c "import Spheral"

    #echo "Spheral ats test..."
    #./spheral-atstest tests/integration.ats

    ## TODO inject debug and non mpi filters when appropriate ...
    #echo "Running spheral-atstest from ${spheral_install_prefix} ..."

    #filter_opt=""

    #if [[ ${spec} == *"~mpi"* ]]; then
    #  fitler_opt="${filter_opt}--filter=\"\'np<2\'\" "
    #  echo "MPI is disabled in the spec running tests with ${filter_opt}"
    #fi

    #spheral-atstest ${spheral_install_prefix}/tests/integration.ats ${filter_opt}

fi
date

