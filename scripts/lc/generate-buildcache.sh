trap 'echo "# $BASH_COMMAND"' DEBUG

### Create a tar file containing:
# dev-pkg/
#   *cloned spheral repo.*
#   resources/
#     pip_cache/
#     mirror/
#     build_cache/
#     bootstrap/
#       metadata/
#       bootstrap_cache/

###############################################################################
###############################################################################

# The expected spack package name for what we are packing up.
SPACK_PKG_NAME=${SPACK_PKG_NAME:-spheral}

# What spec are we targetting.
SPEC=${SPEC:-gcc@10.3.1}

# What is the version of spehral.
SPHERAL_REV_STR=${SPHERAL_REV_STR:-undefined}

# Where will we be staging the package as it is being compiled together.
STAGE_DIR=${STAGE_DIR:-$PWD/../$SYS_TYPE/spheral-cache}

# Where is the spack upstream located.
UPSTREAM_DIR=${UPSTREAM_DIR:-/usr/WS2/sduser/Spheral/spack_upstream/0.22}

# Where does the spheral pip_cache dir live.
SPHERAL_PIP_CACHE_DIR=${SPHERAL_PIP_CACHE_DIR:-~/.cache/spheral_pip}

# What is the local script path.
SCRIPT_DIR=${SCRIPT_DIR:-'scripts'}

# CI_PROJECT_DIR Assuming this is run under Gitlab CI is where a cloned
# version of Spehral lives with no build directories in the dir tree.
CI_PROJECT_DIR=${CI_PROJECT_DIR:-$PWD}

###############################################################################
###############################################################################

# DEV_PKG_NAME defines the title of the spheral install. This is a combination of
# the system type and the current spheral version string.
DEV_PKG_NAME=$SYS_TYPE-spheral-dev-pkg-$SPHERAL_REV_STR

# Full path of what the package directory will look like as we compiler the dev pkg.
DEV_PKG_DIR=$STAGE_DIR/$DEV_PKG_NAME

# Full Spack spec.
SPHERAL_SPEC=$SPACK_PKG_NAME@develop$SPEC

# RESOURCE_DIR is a directory created internally to maintain spack & pip
# resources required for building and running Spheral
RESOURCE_DIR=$DEV_PKG_DIR/resources

# Print for sanity check.
echo $SPHERAL_SPEC
echo $RESOURCE_DIR
echo $STAGE_DIR
echo $DEV_PKG_DIR
echo $SPHERAL_PIP_CACHE_DIR

# Clear the stage directory, create resource dir and copy the Spheral repo into
# the DEV_PKG_DIR.
rm -rf $STAGE_DIR
mkdir -p $RESOURCE_DIR && cp -a $CI_PROJECT_DIR/. $DEV_PKG_DIR

# Copy the SPHERAL_PIP_CACHE_DIR into resource.
mkdir -p $RESOURCE_DIR/pip_cache
cp -a $SPHERAL_PIP_CACHE_DIR/. $RESOURCE_DIR/pip_cache

# Initialize the upstream spack repo.
./$SCRIPT_DIR/devtools/tpl-manager.py --init-only --spack-dir=$UPSTREAM_DIR
source $UPSTREAM_DIR/spack/share/spack/setup-env.sh

# Delete any semblance of a spack env in the STAGE_DIR.
spack env rm -y -f $STAGE_DIR

# Create a spack env in STAGE_DIR and activate it.
spack env create -d $STAGE_DIR
spack env activate $STAGE_DIR

# Concretize our targetted SPHERAL_SPEC.
spack add $SPHERAL_SPEC
spack concretize -f --fresh --deprecated

# Create a mirror of all tpl specs in our environment
# (should only be our deps for SPHERAL_SPEC in the env).
spack mirror create -a -d $RESOURCE_DIR/mirror --exclude-specs "llnlspheral spheral"

# Use spack to list all specs in the mirror and push them to the buildcache.
spack buildcache push -auf $RESOURCE_DIR/mirror $(spack find --format /{hash})

# Mirror bootstrap packages needed to start a spack instance on an airgapped system.
spack bootstrap mirror --binary-packages $RESOURCE_DIR

# Tar up everything in the STAGE_DIR.
tar -czf $DEV_PKG_DIR.tar.gz -C $STAGE_DIR $DEV_PKG_NAME
