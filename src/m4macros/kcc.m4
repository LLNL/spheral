AC_DEFUN([SETUP_KCC], [
  if test -e /usr/local/kcc.base/kcc4.0f/KCC_BASE/bin/KCC; then
    echo "Setting up helpers/myKCC"
    cp /usr/local/kcc.base/kcc4.0f/KCC_BASE/bin/KCC helpers/myKCC
    patch helpers/myKCC helpers/KCC.diff
    rm -f helpers/myKCC.orig
    echo "Setting up helpers/mympKCC"
    rm -f helpers/mympKCC
    CWD=`pwd`
    helpers/processmympKCC.py $CWD
    chmod +x helpers/mympKCC
  else
    echo "Can't find KCC to patch, helpers/myKCC not created."
  fi
])

