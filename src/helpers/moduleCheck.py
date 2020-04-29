from modulefinder import ModuleFinder
import sys, os
import filecmp

mod_name = sys.argv[1]
mod_file = sys.argv[2]

finder = ModuleFinder()
finder.run_script(mod_file)

current_stamp_name = mod_name + "_stamp.cmake"
tmp_stamp_name = current_stamp_name + ".tmp"

newF = open(tmp_stamp_name, "w")
newF.write("set("+mod_name+"_DEPENDS \n")

for name, mod in finder.modules.iteritems():
  if (mod.__file__):
    if not ("lib/python2.7" in mod.__file__):
      newF.write(mod.__file__)
      newF.write('\n')

newF.write(")\n")
newF.close()

if (os.path.isfile(current_stamp_name)):
  if (not filecmp.cmp(current_stamp_name, tmp_stamp_name)):
    os.rename(tmp_stamp_name, current_stamp_name)
  else:
    os.remove(tmp_stamp_name)
else:
  os.rename(tmp_stamp_name, current_stamp_name)
