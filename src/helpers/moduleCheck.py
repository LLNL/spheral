from modulefinder import ModuleFinder
import sys

finder = ModuleFinder()
finder.run_script(sys.argv[1])

out = ""
for name, mod in finder.modules.iteritems():
  if (mod.__file__):
    if not ("lib/python2.7" in mod.__file__):
      out = out + mod.__file__ + ";"
print out
