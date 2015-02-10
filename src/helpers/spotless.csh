#!/usr/bin/env csh
foreach dir (`find . -type d | grep -v 'ti_files' | grep -v 'CVS' | grep -v 'tests'`)
pushd $dir
make spotless
popd
end
