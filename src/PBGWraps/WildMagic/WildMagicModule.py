from pybindgen import *
from PBGutils import *

double = Parameter.new("double", "double")

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class WildMagic:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod):

        # Includes.
        mod.add_include('"WildMagic/WildMagicTypes.hh"')
    
        # Namespace.
        Wm5 = mod.add_cpp_namespace("Wm5")

        # Expose types.
        self.WMVector2d = addObject(Wm5, "WMVector2d")
        self.WMVector3d = addObject(Wm5, "WMVector3d")

        self.WMBox2d = addObject(Wm5, "WMBox2d")
        self.WMBox3d = addObject(Wm5, "WMBox3d")

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.addWMVector2dMethods()
        self.addWMVector3dMethods()

        self.addWMBox2dMethods()
        self.addWMBox3dMethods()

        return

    #---------------------------------------------------------------------------
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["Wm5"]

    #-------------------------------------------------------------------------------
    # Helper method for wrapping WMVector2d
    #-------------------------------------------------------------------------------
    def addWMVector2dMethods(self):
    
        x = self.WMVector2d

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("double", "x"), param("double", "y")])
        
        # Attributes.
        x.add_instance_attribute("X", "double", getter="X", is_const=True)
        x.add_instance_attribute("Y", "double", getter="Y", is_const=True)

        # Operators.
        x.add_unary_numeric_operator("-")
    
        x.add_binary_numeric_operator("+")
        x.add_binary_numeric_operator("-")
    
        x.add_inplace_numeric_operator("+=")
        x.add_inplace_numeric_operator("-=")
    
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")
    
        x.add_binary_numeric_operator("*", left_cppclass = double)
    
        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
    
        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")

        # String representation.
        x.add_function_as_method("printReprWMVector2d", "std::string", 
                                 [param("WMVector2d", "self")],
                                 custom_name = "__repr__")

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping WMVector3d
    #-------------------------------------------------------------------------------
    def addWMVector3dMethods(self):
    
        x = self.WMVector3d

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("double", "x"), param("double", "y"), param("double", "z")])
        
        # Attributes.
        x.add_instance_attribute("X", "double", getter="X", is_const=True)
        x.add_instance_attribute("Y", "double", getter="Y", is_const=True)

        # Operators.
        x.add_unary_numeric_operator("-")
    
        x.add_binary_numeric_operator("+")
        x.add_binary_numeric_operator("-")
    
        x.add_inplace_numeric_operator("+=")
        x.add_inplace_numeric_operator("-=")
    
        x.add_binary_numeric_operator("*", right = "double")
        x.add_binary_numeric_operator("/", right = "double")
    
        x.add_binary_numeric_operator("*", left_cppclass = double)
    
        x.add_inplace_numeric_operator("*=", right = "double")
        x.add_inplace_numeric_operator("/=", right = "double")
    
        # Comparisons.
        x.add_binary_comparison_operator("==")
        x.add_binary_comparison_operator("!=")
        x.add_binary_comparison_operator("<")
        x.add_binary_comparison_operator(">")
        x.add_binary_comparison_operator("<=")
        x.add_binary_comparison_operator(">=")

        # String representation.
        x.add_function_as_method("printReprWMVector3d", "std::string", 
                                 [param("WMVector3d", "self")],
                                 custom_name = "__repr__")

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping WMBox2d
    #-------------------------------------------------------------------------------
    def addWMBox2dMethods(self):
    
        x = self.WMBox2d

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("WMVector2d", "center"),
                           param("WMVector2d", "axis0"),
                           param("WMVector2d", "axis1"),
                           param("double", "extent0"),
                           param("double", "extent1")])

        # Attributes.
        x.add_instance_attribute("Center", "WMVector2d")
        
        # String representation.
        x.add_function_as_method("printReprWMBox2d", "std::string", 
                                 [param("WMBox2d", "self")],
                                 custom_name = "__repr__")

        return

    #-------------------------------------------------------------------------------
    # Helper method for wrapping WMBox3d
    #-------------------------------------------------------------------------------
    def addWMBox3dMethods(self):
    
        x = self.WMBox3d

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("WMVector3d", "center"),
                           param("WMVector3d", "axis0"),
                           param("WMVector3d", "axis1"),
                           param("WMVector3d", "axis2"),
                           param("double", "extent0"),
                           param("double", "extent1"),
                           param("double", "extent2")])

        # Attributes.
        x.add_instance_attribute("Center", "WMVector3d")
        
        # String representation.
        x.add_function_as_method("printReprWMBox3d", "std::string", 
                                 [param("WMBox3d", "self")],
                                 custom_name = "__repr__")

        return

