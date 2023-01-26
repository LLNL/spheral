#-------------------------------------------------------------------------------
# For the given script, find the files it depends on from imports
#-------------------------------------------------------------------------------
import os, sys, modulefinder

# Parse input for the file we're inspecting and any paths we want to exclude
if len(sys.argv) == 2:
    modfile, exclude_paths = sys.argv[1], []
else:
    modfile, exclude_paths = sys.argv[1], sys.argv[2:]

# Define a function we can call recursively for adding up the files
result = []
def addModFile(mfile):
    mfile = mfile.replace(".pyc", ".py")
    if (mfile and
        (os.path.splitext(mfile)[1] == ".py") and
        (not mfile in result) and
        (not max([False] + [mfile.startswith(x) for x in exclude_paths]))):
        result.append(mfile)
        f = modulefinder.ModuleFinder()
        f.run_script(mfile)
        for othername, othermod in f.modules.items():
            if othermod.__file__ != mfile and othermod.__file__ and not max([False] + [othermod.__file__.startswith(x) for x in exclude_paths]):
                addModFile(othermod.__file__)

# Do the deed
addModFile(modfile)

# Output in a make dependency friendly manner
for x in result[1:]:   # Skip the file that started this all off
    print(" %s \\" % x)
