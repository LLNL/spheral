#!/usr/bin/env bash

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
    if [[ "${option}" != "--tpl-only" ]]
    then
        #prefix="/dev/shm/${hostname}"
        prefix="shm/${hostname}"
        if [[ -z ${job_unique_id} ]]; then
          job_unique_id=manual_job_$(date +%s)
          while [[ -d ${prefix}/${job_unique_id} ]] ; do
              sleep 1
              job_unique_id=manual_job_$(date +%s)
          done
        fi

        prefix="${prefix}/${job_unique_id}"
        #prefix="/usr/WS2/davis291/SPHERAL/spheral-pr/shm/rzgenie47/manual_job_1637216151"
        mkdir -p ${prefix}
        prefix_opt="--prefix=${prefix}"

        upstream_opt="--upstream=/usr/WS2/davis291/SPHERAL/lc_uberenv_tpl2/spack/opt/spack"
    fi

    python3 scripts/uberenv/uberenv.py --spec="${spec}" --reuse=True ${upstream_opt} ${prefix_opt}

    echo "Activating spack instance."
    . ${prefix}/spack/share/spack/setup-env.sh

    spheral_spack_hash=`spack find -L spheral | grep spheral | cut -d ' ' -f1`
    spheral_install_prefix=`spack find -p spheral | grep spheral | cut -d ' ' -f3`

    echo "Loading spheral ${spheral_spack_hash} ..."
    spack load /${spheral_spack_hash}

    # TODO inject debug and non mpi filters when appropriate ...
    echo "Running spheral-atstest from ${spheral_install_prefix} ..."

    filter_opt=""

    if [[ ${spec} == *"~mpi"* ]]; then
      fitler_opt="${filter_opt}--filter=\"\'np<2\'\" "
      echo "MPI is disabled in the spec running tests with ${filter_opt}"
    fi

    spheral-atstest ${spheral_install_prefix}/tests/integration.ats ${filter_opt}

fi
date

