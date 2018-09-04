from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class FieldOperations:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):
        self.dims = dims
        mod.add_include('"%s/FieldOperationsTypes.hh"' % srcdir)
        self.space = mod.add_cpp_namespace("Spheral")
        return

    #---------------------------------------------------------------------------
    # Generate bindings.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        for dim in self.dims:
            self.generateDimBindings(dim)
        return

    #---------------------------------------------------------------------------
    # Add the types per dimension.
    #---------------------------------------------------------------------------
    def generateDimBindings(self, ndim):
        space = self.space

        dim = "Spheral::Dim<%i>" % ndim
        vector = "Spheral::Vector%id" % ndim
        tensor = "Spheral::Tensor%id" % ndim
        symtensor = "Spheral::SymTensor%id" % ndim
        intfieldlist = "Spheral::IntFieldList%id" % ndim
        scalarfieldlist = "Spheral::ScalarFieldList%id" % ndim
        vectorfieldlist = "Spheral::VectorFieldList%id" % ndim
        tensorfieldlist = "Spheral::TensorFieldList%id" % ndim
        symtensorfieldlist = "Spheral::SymTensorFieldList%id" % ndim
        fieldlistset = "Spheral::FieldListSet%id" % ndim
        tablekernel = "Spheral::TableKernel%id" % ndim
        vector_of_vector_of_Vector = "vector_of_vector_of_Vector%id" % ndim
        vector_of_vector_of_Tensor = "vector_of_vector_of_Tensor%id" % ndim
        vector_of_vector_of_SymTensor = "vector_of_vector_of_SymTensor%id" % ndim
        vector_of_Boundary = "vector_of_Boundary%id" % ndim

        # smoothFields
        space.add_function("smoothFields", scalarfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                             constrefparam(vectorfieldlist, "position"),
                                                             constrefparam(scalarfieldlist, "weight"),
                                                             constrefparam(scalarfieldlist, "mass"),
                                                             constrefparam(scalarfieldlist, "rho"),
                                                             constrefparam(symtensorfieldlist, "Hfield"),
                                                             constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, "double"],
                           custom_name = "smoothScalarFields%id" % ndim)
        space.add_function("smoothFields", vectorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                             constrefparam(vectorfieldlist, "position"),
                                                             constrefparam(scalarfieldlist, "weight"),
                                                             constrefparam(scalarfieldlist, "mass"),
                                                             constrefparam(scalarfieldlist, "rho"),
                                                             constrefparam(symtensorfieldlist, "Hfield"),
                                                             constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "smoothVectorFields%id" % ndim)
        space.add_function("smoothFields", tensorfieldlist, [constrefparam(tensorfieldlist, "fieldList"),
                                                             constrefparam(vectorfieldlist, "position"),
                                                             constrefparam(scalarfieldlist, "weight"),
                                                             constrefparam(scalarfieldlist, "mass"),
                                                             constrefparam(scalarfieldlist, "rho"),
                                                             constrefparam(symtensorfieldlist, "Hfield"),
                                                             constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, tensor],
                           custom_name = "smoothTensorFields%id" % ndim)
        space.add_function("smoothFields", symtensorfieldlist, [constrefparam(symtensorfieldlist, "fieldList"),
                                                                constrefparam(vectorfieldlist, "position"),
                                                                constrefparam(scalarfieldlist, "weight"),
                                                                constrefparam(scalarfieldlist, "mass"),
                                                                constrefparam(scalarfieldlist, "rho"),
                                                                constrefparam(symtensorfieldlist, "Hfield"),
                                                                constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, symtensor],
                           custom_name = "smoothSymTensorFields%id" % ndim)

        # smoothFieldsMash
        space.add_function("smoothFieldsMash", scalarfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, "double"],
                           custom_name = "smoothScalarFieldsMash%id" % ndim)
        space.add_function("smoothFieldsMash", vectorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "smoothVectorFieldsMash%id" % ndim)
        space.add_function("smoothFieldsMash", tensorfieldlist, [constrefparam(tensorfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, tensor],
                           custom_name = "smoothTensorFieldsMash%id" % ndim)
        space.add_function("smoothFieldsMash", symtensorfieldlist, [constrefparam(symtensorfieldlist, "fieldList"),
                                                                    constrefparam(vectorfieldlist, "position"),
                                                                    constrefparam(scalarfieldlist, "weight"),
                                                                    constrefparam(symtensorfieldlist, "Hfield"),
                                                                    constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, symtensor],
                           custom_name = "smoothSymTensorFieldsMash%id" % ndim)

        # sampleFieldsMash
        space.add_function("sampleFieldsMash", scalarfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel"),
                                                                 constrefparam(vectorfieldlist, "samplePosition"),
                                                                 constrefparam(scalarfieldlist, "sampleWeight"),
                                                                 constrefparam(symtensorfieldlist, "sampleHfield")],
                           template_parameters = [dim, "double"],
                           custom_name = "sampleScalarFieldsMash%id" % ndim)
        space.add_function("sampleFieldsMash", vectorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel"),
                                                                 constrefparam(vectorfieldlist, "samplePosition"),
                                                                 constrefparam(scalarfieldlist, "sampleWeight"),
                                                                 constrefparam(symtensorfieldlist, "sampleHfield")],
                           template_parameters = [dim, vector],
                           custom_name = "sampleVectorFieldsMash%id" % ndim)
        space.add_function("sampleFieldsMash", tensorfieldlist, [constrefparam(tensorfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel"),
                                                                 constrefparam(vectorfieldlist, "samplePosition"),
                                                                 constrefparam(scalarfieldlist, "sampleWeight"),
                                                                 constrefparam(symtensorfieldlist, "sampleHfield")],
                           template_parameters = [dim, tensor],
                           custom_name = "sampleTensorFieldsMash%id" % ndim)
        space.add_function("sampleFieldsMash", symtensorfieldlist, [constrefparam(symtensorfieldlist, "fieldList"),
                                                                    constrefparam(vectorfieldlist, "position"),
                                                                    constrefparam(scalarfieldlist, "weight"),
                                                                    constrefparam(symtensorfieldlist, "Hfield"),
                                                                    constrefparam(tablekernel, "kernel"),
                                                                    constrefparam(vectorfieldlist, "samplePosition"),
                                                                    constrefparam(scalarfieldlist, "sampleWeight"),
                                                                    constrefparam(symtensorfieldlist, "sampleHfield")],
                           template_parameters = [dim, symtensor],
                           custom_name = "sampleSymTensorFieldsMash%id" % ndim)

        # sampleMultipleFieldsMash
        space.add_function("sampleMultipleFieldsMash", fieldlistset, [constrefparam(fieldlistset, "fieldListSet"),
                                                                      constrefparam(vectorfieldlist, "position"),
                                                                      constrefparam(scalarfieldlist, "weight"),
                                                                      constrefparam(symtensorfieldlist, "Hfield"),
                                                                      constrefparam(tablekernel, "kernel"),
                                                                      constrefparam(vectorfieldlist, "samplePosition"),
                                                                      constrefparam(scalarfieldlist, "sampleWeight"),
                                                                      constrefparam(symtensorfieldlist, "sampleHfield")],
                           template_parameters = [dim],
                           custom_name = "sampleMultipleScalarFieldsMash%id" % ndim)

        # splatFieldsMash
        space.add_function("splatFieldsMash", scalarfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                                constrefparam(vectorfieldlist, "position"),
                                                                constrefparam(scalarfieldlist, "weight"),
                                                                constrefparam(symtensorfieldlist, "Hfield"),
                                                                constrefparam(tablekernel, "kernel"),
                                                                constrefparam(vectorfieldlist, "splatPosition"),
                                                                constrefparam(scalarfieldlist, "splatWeight"),
                                                                constrefparam(symtensorfieldlist, "splatHfield")],
                           template_parameters = [dim, "double"],
                           custom_name = "splatScalarFieldsMash%id" % ndim)
        space.add_function("splatFieldsMash", vectorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                                constrefparam(vectorfieldlist, "position"),
                                                                constrefparam(scalarfieldlist, "weight"),
                                                                constrefparam(symtensorfieldlist, "Hfield"),
                                                                constrefparam(tablekernel, "kernel"),
                                                                constrefparam(vectorfieldlist, "splatPosition"),
                                                                constrefparam(scalarfieldlist, "splatWeight"),
                                                                constrefparam(symtensorfieldlist, "splatHfield")],
                           template_parameters = [dim, vector],
                           custom_name = "splatVectorFieldsMash%id" % ndim)
        space.add_function("splatFieldsMash", tensorfieldlist, [constrefparam(tensorfieldlist, "fieldList"),
                                                                constrefparam(vectorfieldlist, "position"),
                                                                constrefparam(scalarfieldlist, "weight"),
                                                                constrefparam(symtensorfieldlist, "Hfield"),
                                                                constrefparam(tablekernel, "kernel"),
                                                                constrefparam(vectorfieldlist, "splatPosition"),
                                                                constrefparam(scalarfieldlist, "splatWeight"),
                                                                constrefparam(symtensorfieldlist, "splatHfield")],
                           template_parameters = [dim, tensor],
                           custom_name = "splatTensorFieldsMash%id" % ndim)
        space.add_function("splatFieldsMash", symtensorfieldlist, [constrefparam(symtensorfieldlist, "fieldList"),
                                                                   constrefparam(vectorfieldlist, "position"),
                                                                   constrefparam(scalarfieldlist, "weight"),
                                                                   constrefparam(symtensorfieldlist, "Hfield"),
                                                                   constrefparam(tablekernel, "kernel"),
                                                                   constrefparam(vectorfieldlist, "splatPosition"),
                                                                   constrefparam(scalarfieldlist, "splatWeight"),
                                                                   constrefparam(symtensorfieldlist, "splatHfield")],
                           template_parameters = [dim, symtensor],
                           custom_name = "splatSymTensorFieldsMash%id" % ndim)

        # splatMultipleFieldsMash
        space.add_function("splatMultipleFieldsMash", fieldlistset, [constrefparam(fieldlistset, "fieldListSet"),
                                                                     constrefparam(vectorfieldlist, "position"),
                                                                     constrefparam(scalarfieldlist, "weight"),
                                                                     constrefparam(symtensorfieldlist, "Hfield"),
                                                                     constrefparam(tablekernel, "kernel"),
                                                                     constrefparam(vectorfieldlist, "splatPosition"),
                                                                     constrefparam(scalarfieldlist, "splatWeight"),
                                                                     constrefparam(symtensorfieldlist, "splatHfield"),
                                                                     constrefparam(vector_of_Boundary, "boundaries")],
                           template_parameters = [dim],
                           custom_name = "splatMultipleFieldsMash%id" % ndim)

        # gradient
        space.add_function("gradient", vectorfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                         constrefparam(vectorfieldlist, "position"),
                                                         constrefparam(scalarfieldlist, "weight"),
                                                         constrefparam(scalarfieldlist, "mass"),
                                                         constrefparam(scalarfieldlist, "rho"),
                                                         constrefparam(symtensorfieldlist, "Hfield"),
                                                         constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, "double"],
                           custom_name = "gradient")
        space.add_function("gradient", tensorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                         constrefparam(vectorfieldlist, "position"),
                                                         constrefparam(scalarfieldlist, "weight"),
                                                         constrefparam(scalarfieldlist, "mass"),
                                                         constrefparam(scalarfieldlist, "rho"),
                                                         constrefparam(symtensorfieldlist, "Hfield"),
                                                         constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "gradient")

        # divergence
        space.add_function("divergence", scalarfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                           constrefparam(vectorfieldlist, "position"),
                                                           constrefparam(scalarfieldlist, "weight"),
                                                           constrefparam(scalarfieldlist, "mass"),
                                                           constrefparam(scalarfieldlist, "rho"),
                                                           constrefparam(symtensorfieldlist, "Hfield"),
                                                           constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "divergence")
        space.add_function("divergence", vectorfieldlist, [constrefparam(tensorfieldlist, "fieldList"),
                                                           constrefparam(vectorfieldlist, "position"),
                                                           constrefparam(scalarfieldlist, "weight"),
                                                           constrefparam(scalarfieldlist, "mass"),
                                                           constrefparam(scalarfieldlist, "rho"),
                                                           constrefparam(symtensorfieldlist, "Hfield"),
                                                           constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, tensor],
                           custom_name = "divergence")
        space.add_function("divergence", vectorfieldlist, [constrefparam(symtensorfieldlist, "fieldList"),
                                                           constrefparam(vectorfieldlist, "position"),
                                                           constrefparam(scalarfieldlist, "weight"),
                                                           constrefparam(scalarfieldlist, "mass"),
                                                           constrefparam(scalarfieldlist, "rho"),
                                                           constrefparam(symtensorfieldlist, "Hfield"),
                                                           constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, symtensor],
                           custom_name = "divergence")

        # gradientMash
        space.add_function("gradientMash", vectorfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                             constrefparam(vectorfieldlist, "position"),
                                                             constrefparam(scalarfieldlist, "weight"),
                                                             constrefparam(symtensorfieldlist, "Hfield"),
                                                             constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, "double"],
                           custom_name = "gradientMash")
        space.add_function("gradientMash", tensorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                             constrefparam(vectorfieldlist, "position"),
                                                             constrefparam(scalarfieldlist, "weight"),
                                                             constrefparam(symtensorfieldlist, "Hfield"),
                                                             constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "gradientMash")

        # gradientPairWise
        space.add_function("gradientPairWise", vectorfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(scalarfieldlist, "mass"),
                                                                 constrefparam(scalarfieldlist, "rho"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, "double"],
                           custom_name = "gradientPairWise")
        space.add_function("gradientPairWise", tensorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                                 constrefparam(vectorfieldlist, "position"),
                                                                 constrefparam(scalarfieldlist, "weight"),
                                                                 constrefparam(scalarfieldlist, "mass"),
                                                                 constrefparam(scalarfieldlist, "rho"),
                                                                 constrefparam(symtensorfieldlist, "Hfield"),
                                                                 constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "gradientPairWise")

        # sampleMultipleFields2Lattice
        space.add_function("sampleMultipleFields2Lattice", None, [constrefparam(fieldlistset, "fieldListSet"),
                                                                  constrefparam(vectorfieldlist, "position"),
                                                                  constrefparam(scalarfieldlist, "weight"),
                                                                  constrefparam(symtensorfieldlist, "Hfield"),
                                                                  constrefparam(intfieldlist, "mask"),
                                                                  constrefparam(tablekernel, "kernel"),
                                                                  constrefparam(vector, "xmin"),
                                                                  constrefparam(vector, "xmax"),
                                                                  constrefparam("vector_of_int", "nsample"),
                                                                  refparam("vector_of_vector_of_double", "scalarValues"),
                                                                  refparam(vector_of_vector_of_Vector, "vectorValues"),
                                                                  refparam(vector_of_vector_of_Tensor, "tensorValues"),
                                                                  refparam(vector_of_vector_of_SymTensor, "symTensorValues")],
                           template_parameters = [dim],
                           custom_name = "sampleMultipleFields2Lattice%id" % ndim)

        # sampleMultipleFields2LatticeMash
        space.add_function("sampleMultipleFields2LatticeMash", None, [constrefparam(fieldlistset, "fieldListSet"),
                                                                      constrefparam(vectorfieldlist, "position"),
                                                                      constrefparam(scalarfieldlist, "weight"),
                                                                      constrefparam(symtensorfieldlist, "Hfield"),
                                                                      constrefparam(intfieldlist, "mask"),
                                                                      constrefparam(tablekernel, "kernel"),
                                                                      constrefparam(vector, "xmin"),
                                                                      constrefparam(vector, "xmax"),
                                                                      constrefparam("vector_of_int", "nsample"),
                                                                      refparam("vector_of_vector_of_double", "scalarValues"),
                                                                      refparam(vector_of_vector_of_Vector, "vectorValues"),
                                                                      refparam(vector_of_vector_of_Tensor, "tensorValues"),
                                                                      refparam(vector_of_vector_of_SymTensor, "symTensorValues")],
                           template_parameters = [dim],
                           custom_name = "sampleMultipleFields2LatticeMash%id" % ndim)

        # limiter
        space.add_function("limiter", symtensorfieldlist, [constrefparam(scalarfieldlist, "fieldList"),
                                                           constrefparam(vectorfieldlist, "gradient"),
                                                           constrefparam(vectorfieldlist, "position"),
                                                           constrefparam(symtensorfieldlist, "Hfield"),
                                                           constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, "double"],
                           custom_name = "scalarLimiter%id" % ndim)
        space.add_function("limiter", symtensorfieldlist, [constrefparam(vectorfieldlist, "fieldList"),
                                                           constrefparam(tensorfieldlist, "gradient"),
                                                           constrefparam(vectorfieldlist, "position"),
                                                           constrefparam(symtensorfieldlist, "Hfield"),
                                                           constrefparam(tablekernel, "kernel")],
                           template_parameters = [dim, vector],
                           custom_name = "vectorLimiter%id" % ndim)

        # binFieldList2Lattice
        space.add_function("binFieldList2Lattice", None, [constrefparam(scalarfieldlist, "fieldList"),
                                                          constrefparam(vector, "xmin"),
                                                          constrefparam(vector, "xmax"),
                                                          constrefparam("vector_of_unsigned", "nsample")],
                           template_parameters = [dim, "double"],
                           custom_name = "binScalarFieldList2Lattice%id" % ndim)
        space.add_function("binFieldList2Lattice", None, [constrefparam(vectorfieldlist, "fieldList"),
                                                          constrefparam(vector, "xmin"),
                                                          constrefparam(vector, "xmax"),
                                                          constrefparam("vector_of_unsigned", "nsample")],
                           template_parameters = [dim, vector],
                           custom_name = "binVectorFieldList2Lattice%id" % ndim)
        space.add_function("binFieldList2Lattice", None, [constrefparam(tensorfieldlist, "fieldList"),
                                                          constrefparam(vector, "xmin"),
                                                          constrefparam(vector, "xmax"),
                                                          constrefparam("vector_of_unsigned", "nsample")],
                           template_parameters = [dim, tensor],
                           custom_name = "binTensorFieldList2Lattice%id" % ndim)
        space.add_function("binFieldList2Lattice", None, [constrefparam(symtensorfieldlist, "fieldList"),
                                                          constrefparam(vector, "xmin"),
                                                          constrefparam(vector, "xmax"),
                                                          constrefparam("vector_of_unsigned", "nsample")],
                           template_parameters = [dim, symtensor],
                           custom_name = "binSymTensorFieldList2Lattice%id" % ndim)

        # binFieldList2Lattice (with kernel smoothing)
        space.add_function("binFieldList2LatticeWithSmoothing", None, [constrefparam(scalarfieldlist, "fieldList"),
                                                                       constrefparam(tablekernel, "W"),
                                                                       constrefparam(vector, "xmin"),
                                                                       constrefparam(vector, "xmax"),
                                                                       constrefparam("vector_of_unsigned", "nsample")],
                           template_parameters = [dim, "double"],
                           custom_name = "binScalarFieldList2LatticeWithSmoothing%id" % ndim)
        return

