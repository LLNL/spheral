#-------------------------------------------------------------------------------
# For the given script, find the files it depends on from imports
#-------------------------------------------------------------------------------
import os, sys, modulefinder

# Parse input for the file we're inspecting and any paths we want to exclude
if len(sys.argv) == 2:
    modfile, exclude_paths = sys.argv[1], []
else:
    modfile, exclude_paths = sys.argv[1], sys.argv[1:]

#exec("import %s as mod" % os.path.splitext(modfile)[0])
#modnames, thpt = externals(mod) # os.path.splitext(modfile)[0])

# Define a function we can call recursively for adding up the files
result = []
def addModFile(mod):
    mfile = mod.__file__.replace(".pyc", ".py")
    if mfile and (not mfile in result) and (os.path.splitext(mfile)[1] == ".py"):
        ok = not max([False] + [mfile.startswith(x) for x in exclude_paths])
        if ok:
            result.append(mfile)
            f = modulefinder.ModuleFinder()
            f.run_script(mfile)
            for othername, othermod in f.modules.iteritems():
                if othermod.__file__ != mfile and othermod.__file__ and not max([False] + [othermod.__file__.startswith(x) for x in exclude_paths]):
                    addModFile(othermod)

# Do the deed
exec("import %s as mod" % os.path.splitext(modfile)[0])
addModFile(mod)

for x in result:
    print x
