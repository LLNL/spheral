# -*- mode: python -*-
#-------------------------------------------------------------------------------
# Provide a method of getting the instantiated dimensions from python.
#-------------------------------------------------------------------------------
def spheralDimensions():
    spheralDimensions.PYB11ignore = True     # Screen from PYB11
    dims_string = "@DIMS@"
    dims = [int(x) for x in dims_string.split()]
    return dims

#-------------------------------------------------------------------------------
# Return a dictionary of useful Spheral Dimensional types.
#-------------------------------------------------------------------------------
def dimDictionary(ndim):
    dimDictionary.PYB11ignore = True     # Screen from PYB11
    return {"Dimension"      : "Dim<%i>" % ndim,
            "Scalar"         : "Dim<%i>::Scalar" % ndim,
            "Vector"         : "Dim<%i>::Scalar" % ndim,
            "Tensor"         : "Dim<%i>::Scalar" % ndim,
            "SymTensor"      : "Dim<%i>::Scalar" % ndim,
            "ScalarField"    : "Field<Dim<%i>, Dim<%i>::Scalar>" % (ndim, ndim),
            "VectorField"    : "Field<Dim<%i>, Dim<%i>::Vector>" % (ndim, ndim),
            "TensorField"    : "Field<Dim<%i>, Dim<%i>::Tensor>" % (ndim, ndim),
            "SymTensorField" : "Field<Dim<%i>, Dim<%i>::SymTensor>" % (ndim, ndim),
            }
dimDictionary.PYB11ignore = True     # Screen from PYB11
