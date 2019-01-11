from pybindgen import *

from PBGutils import *
from ref_return_value import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class FractalGravity:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/FractalGravity/FractalGravity.hh"' % topsrcdir)
    
        # Namespace.
        space = mod.add_cpp_namespace("Spheral")
        genericbodyforce3d = findObject(space, "GenericBodyForce3d")

        # Expose types.
        self.FractalGravity = addObject(space, "FractalGravity", allow_subclassing=True, parent=genericbodyforce3d)

        return

    #---------------------------------------------------------------------------
    # Generate bindings for all objects.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):
        self.generateFractalGravityBindings(self.FractalGravity)
        return

    #---------------------------------------------------------------------------
    # FractalGravity bindings.
    #---------------------------------------------------------------------------
    def generateFractalGravityBindings(self, x):

        # Object names.
        me = "Spheral::FractalGravity"
        state = "Spheral::State3d"
        derivatives = "Spheral::StateDerivatives3d"
        database = "Spheral::DataBase3d"
        scalarfieldlist = "Spheral::ScalarFieldList3d"

        # Constructors.
        x.add_constructor([param("double", "G"),
                           constrefparam("Vector3d", "xmin"),
                           constrefparam("Vector3d", "xmax"),
                           param("bool", "periodic"),
                           param("unsigned int", "ngrid"),
                           param("unsigned int", "nlevelmax"),
                           param("unsigned int", "minHighParticles"),
                           param("unsigned int", "padding"),
                           param("double", "maxDeltaVelocity")])

        # Methods.
        x.add_method("evaluateDerivatives", None, [param("const double", "time"),
                                                   param("const double", "dt"),
                                                   constrefparam(database, "dataBase"),
                                                   constrefparam(state, "state"),
                                                   refparam(derivatives, "derivatives")],
                     is_const=True, is_virtual=True)
        x.add_method("dt", "pair_double_string", [constrefparam(database, "dataBase"),
                                                  constrefparam(state, "state"),
                                                  constrefparam(derivatives, "derivatives"),
                                                  param("double", "time")],
                     is_const=True, is_virtual=True)
        x.add_method("initializeProblemStartup", None, [refparam(database, "dataBase")], is_virtual=True)
        const_ref_return_value(x, me, "%s::potential" % me, scalarfieldlist, [], "potential")

        # Attributes.
        x.add_instance_attribute("G", "double", getter="G", is_const=True)
        x.add_instance_attribute("xmin", "Vector3d", getter="xmin", is_const=True)
        x.add_instance_attribute("xmax", "Vector3d", getter="xmax", is_const=True)
        x.add_instance_attribute("periodic", "bool", getter="periodic", is_const=True)
        x.add_instance_attribute("ngrid", "unsigned int", getter="ngrid", is_const=True)
        x.add_instance_attribute("nlevelmax", "unsigned int", getter="nlevelmax", is_const=True)
        x.add_instance_attribute("minHighParticles", "unsigned int", getter="minHighParticles", is_const=True)
        x.add_instance_attribute("padding", "unsigned int", getter="padding", is_const=True)

        return
