#-------------------------------------------------------------------------------
# For the given script, find the files it depends on from imports
#-------------------------------------------------------------------------------
import os, sys, modulefinder
#from pydeps.pydeps import externals

print sys.argv

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
    print " --> ", mod, mfile
    if mfile:
        ok = not max([False] + [mfile.startswith(x) for x in exclude_paths])
        if ok:
            result.append(mfile)
            f = modulefinder.ModuleFinder()
            f.run_script(mfile)
            for othername, othermod in f.modules.iteritems():
                #print "Check : ", othermod.__file__
                if othermod.__file__ and not max([False] + [othermod.__file__.startswith(x) for x in exclude_paths]):
                    addModFile(othermod)

# Do the deed
print " **> ", os.path.splitext(modfile)[0]
exec("import %s as mod" % os.path.splitext(modfile)[0])
addModFile(mod)

print result
