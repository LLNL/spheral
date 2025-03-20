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

# Where will we be staging the package as it is being compiled together.
INSTALL_DIR=${INSTALL_DIR:-$PWD/../$SYS_TYPE/spheral-cache}

# Where does the spheral pip_cache dir live.
SPHERAL_PIP_CACHE_DIR=${SPHERAL_PIP_CACHE_DIR:-~/.cache/spheral_pip}

# What is the local script path.
SCRIPT_DIR=${SCRIPT_DIR:-'scripts'}

# CI_PROJECT_DIR Assuming this is run under Gitlab CI is where a cloned
# version of Spehral lives with no build directories in the dir tree.
CI_PROJECT_DIR=${CI_PROJECT_DIR:-$PWD}

# DEV_PKG_NAME defines the title of the spheral install. This is a combination of
# the system type and the current spheral version string.
DEV_PKG_NAME=${DEV_PKG_NAME:-$SYS_TYPE-spheral-dev-pkg-undefined}

###############################################################################
###############################################################################

# Full path of what the package directory will look like as we compile the dev pkg.
DEV_PKG_DIR=$INSTALL_DIR/$DEV_PKG_NAME

# RESOURCE_DIR is a directory created internally to maintain spack & pip
# resources required for building and running Spheral
RESOURCE_DIR=$DEV_PKG_DIR/resources

# Print for sanity check.
echo $RESOURCE_DIR
echo $INSTALL_DIR
echo $DEV_PKG_DIR
echo $SPHERAL_PIP_CACHE_DIR

# Clear the stage directory, create resource dir and copy the Spheral repo into
# the DEV_PKG_DIR.
rm -rf $INSTALL_DIR
mkdir -p $RESOURCE_DIR && cp -a $CI_PROJECT_DIR/. $DEV_PKG_DIR

# Copy the SPHERAL_PIP_CACHE_DIR into resource.
mkdir -p $RESOURCE_DIR/pip_cache
cp -a $SPHERAL_PIP_CACHE_DIR/. $RESOURCE_DIR/pip_cache

# tpl-manager --dev-pkg does the following:
# Creates a local Spack repo
# Activates and concretizes the dev_pkg Spheral Spack environment
# Installs the Spheral dependencies for all specs
./$SCRIPT_DIR/devtools/tpl-manager.py --dev-pkg

# Source Spack for the current terminal
source ../spheral-spack-tpls/spack/share/spack/setup-env.sh

# Activate our dev spack environment
spack env activate ./scripts/spack/environments/dev_pkg

# Create a mirror of all tpl specs in our environment
# (should only be our deps for a single spec in the env).
spack mirror create -a -d $RESOURCE_DIR/mirror --exclude-specs "llnlspheral spheral"

# Use spack to list all specs in the mirror and push them to the buildcache.
spack buildcache push -auf $RESOURCE_DIR/mirror $(spack find --format /{hash})

# Mirror bootstrap packages needed to start a spack instance on an airgapped system.
spack bootstrap mirror --binary-packages $RESOURCE_DIR

# Tar up everything in the INSTALL_DIR.
tar -czf $DEV_PKG_DIR.tar.gz -C $INSTALL_DIR $DEV_PKG_NAME
