trap 'echo "# $BASH_COMMAND"' DEBUG

SPACK_PKG_NAME=${SPACK_PKG_NAME:-'spheral'}
SPACK_URL=${SPACK_URL:-'https://github.com/spack/spack'}
BUILD_ALLOC=${BUILD_ALLOC}
SCRIPT_DIR=${SCRIPT_DIR:-'scripts'}

if [[ -z "${SPEC}" ]]; then
  echo "SPEC var must be set."
  exit 1
fi

if [[ -z "${INSTALL_DIR}" ]]; then
  echo "INSTALL_DIR var must be set."
  exit 1
fi

echo $SPACK_PKG_NAME
echo $SPEC
echo $SPACK_URL
echo $INSTALL_DIR
echo $SCRIPT_DIR
echo $BUILD_ALLOC

rm -rf $INSTALL_DIR
mkdir -p $INSTALL_DIR

./$SCRIPT_DIR/devtools/tpl-manager.py --spack-url $SPACK_URL --init-only --spec=none --no-upstream --spheral-spack-dir $INSTALL_DIR/spheral-spack-tpls

source $INSTALL_DIR/spheral-spack-tpls/spack/share/spack/setup-env.sh
spack bootstrap add --trust local-sources $PWD/resources/metadata/sources
spack bootstrap add --trust local-binaries $PWD/resources/metadata/binaries
spack mirror rm spheral-mirror
spack mirror rm spheral-cache
spack mirror add --unsigned spheral-mirror $PWD/resources/mirror
spack mirror add --unsigned spheral-cache $PWD/resources
spack buildcache update-index $PWD/resources/mirror

$BUILD_ALLOC spack install --fresh --deprecated --no-check-signature --only dependencies $SPACK_PKG_NAME@develop%$SPEC

$BUILD_ALLOC ./$SCRIPT_DIR/devtools/tpl-manager.py --spack-url $SPACK_URL --no-upstream --spheral-spack-dir $INSTALL_DIR/spheral-spack-tpls --spec $SPEC
#$BUILD_ALLOC spack install --fresh --deprecated --no-check-signature -u initconfig "$SPACK_PKG_NAME@develop%$SPEC dev_path=$PWD"

HOST_CONFIG_FILE=$(ls -t | grep -E "*\.cmake" | head -1)
$BUILD_ALLOC ./$SCRIPT_DIR/devtools/host-config-build.py --host-config $HOST_CONFIG_FILE -i $INSTALL_DIR --build --no-clean





