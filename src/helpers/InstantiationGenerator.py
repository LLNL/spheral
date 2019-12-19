#-------------------------------------------------------------------------------
# IntantiationGenerator
# 
# A python script to automatically generate Spheral++ instantion files to be 
# compiled.  Assumed arguments:
#    infile - the file to be read, defining "text"
#   outfile - the file to be written out
#      ndim - an integer value for the dimensionality being generated (1,2,3)
#-------------------------------------------------------------------------------
import sys

assert len(sys.argv) == 4
infile = sys.argv[1]
outfile = sys.argv[2]
ndim = sys.argv[3]
idim = int(ndim)

dictionary = {"ndim" : ndim}

## Python 2
#
#execfile(infile)
#
## Python 3
exec(open(infile).read())

f = open(outfile, "w")
f.write(text % dictionary)
f.close()
