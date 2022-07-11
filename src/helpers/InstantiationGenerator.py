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

dictionary = {"ndim"      : ndim,
              "Dim"       : "Dim<%s>" % ndim,
              "Scalar"    : "Dim<%s>::Scalar" % ndim,
              "Vector"    : "Dim<%s>::Vector" % ndim,
              "Tensor"    : "Dim<%s>::Tensor" % ndim,
              "SymTensor" : "Dim<%s>::SymTensor" % ndim,
}

# Read the input file to get the definition of the string "text", which we use to generate the explicit instantiation .cc file
exec(open(infile).read())
with open(outfile, "w") as f:
    f.write(text % dictionary)
