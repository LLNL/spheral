#!/usr/bin/env bash

###############################################################################
# Copyright (c) 2016-22, Lawrence Livermore National Security, LLC
# and Spheral project contributors. See the Spheral/LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

TAR_CMD=`which gtar`
VERSION=`git describe --tags --always`

git archive --prefix=Spheral-${VERSION}/ -o Spheral-${VERSION}.tar HEAD 2> /dev/null

echo "Running git archive submodules..."
echo "Note : If running on OSX you may need to \`brew install gnu-tar\` and change \`TAR_CMD\` to \`which gtar\`..."

p=`pwd` && (echo .; git submodule foreach --recursive) | while read entering path; do
    temp="${path%\'}";
    temp="${temp#\'}";
    path=$temp;
    [ "$path" = "" ] && continue;
    (cd $path && git archive --prefix=Spheral-${VERSION}/$path/ HEAD > $p/tmp.tar && ${TAR_CMD} --concatenate --file=$p/Spheral-${VERSION}.tar $p/tmp.tar && rm $p/tmp.tar);
done

gzip Spheral-${VERSION}.tar
