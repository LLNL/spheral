from PYB11Generator import *

#-------------------------------------------------------------------------------
# GeometryRegistrar
#-------------------------------------------------------------------------------
@PYB11singleton
class GeometryRegistrar:

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11returnpolicy("reference")
    def instance(self):
        return "GeometryRegistrar&"

    # The coordinate system
    @PYB11static
    @PYB11pycppname("coords")
    @PYB11returnpolicy("reference")
    #@PYB11ignore
    def get_coords(self):
        return "CoordinateType"

    @PYB11static
    @PYB11pycppname("coords")
    #@PYB11ignore
    def set_coords(self,
                   x = "const CoordinateType"):
        return "void"

    #coords = property(get_coords, set_coords, doc="The coordinate system")
