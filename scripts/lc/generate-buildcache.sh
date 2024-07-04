trap 'echo "# $BASH_COMMAND"' DEBUG

SPEC=${SPEC:-gcc@10.3.1}
SPACK_PKG_NAME=${SPACK_PKG_NAME:-spheral}

SPHERAL_SPEC=$SPACK_PKG_NAME@develop%$SPEC
echo $SPHERAL_SPEC

UPSTREAM_DIR=${UPSTREAM_DIR:-/usr/WS2/sduser/Spheral/spack_upstream/0.22}
SCRIPT_DIR=${SCRIPT_DIR:-'scripts'}
SPHERAL_REV_STRING=${SPHERAL_REV_STRING:-undefined}
INSTALL_DIR=${INSTALL_DIR:-$PWD/../$SYS_TYPE/spheral-cache}
DEV_PKG_NAME=$SYS_TYPE-spheral-dev-pkg-$SPHERAL_REV_STRING
DEV_PKG_DIR=$INSTALL_DIR/$DEV_PKG_NAME

CI_PROJECT_DIR=${CI_PROJECT_DIR:-$PWD}

RESOURCE_DIR=$DEV_PKG_DIR/resources
echo $RESOURCE_DIR

echo $INSTALL_DIR
#echo $INSTALL_DIR &> install-dir.txt

echo $DEV_PKG_DIR
#echo $DEV_PKG_NAME &> dev-pkg-name.txt

rm -rf $INSTALL_DIR
mkdir -p $RESOURCE_DIR && cp -a $CI_PROJECT_DIR/. $DEV_PKG_DIR

./$SCRIPT_DIR/devtools/tpl-manager.py --init-only --spheral-spack-dir=$UPSTREAM_DIR --spec=none
source $UPSTREAM_DIR/spack/share/spack/setup-env.sh

spack env rm -y -f $INSTALL_DIR
spack env create -d $INSTALL_DIR
spack env activate $INSTALL_DIR
spack add $SPHERAL_SPEC
spack concretize --fresh -f

spack mirror create -a -d $RESOURCE_DIR/mirror --exclude-specs "llnlspheral spheral"
spack mirror set --unsigned $RESOURCE_DIR/mirror

spack buildcache push -auf $RESOURCE_DIR/mirror $(spack find --format /{hash})

spack bootstrap mirror --binary-packages $RESOURCE_DIR

tar -czvf $DEV_PKG_DIR.tar.gz -C $INSTALL_DIR $DEV_PKG_NAME

