from pybindgen import *
from PBGutils import *

#-------------------------------------------------------------------------------
# ref_return_value
#-------------------------------------------------------------------------------
def ref_return_value(self,                 # The class we're adding to
                     selfType,             # The name of the class we're adding to
                     method,               # The accessor method returning a reference
                     returnType,           # The return type (without the &)
                     parameters = [],      # optional : list or tuple of paramters for the method.
                     custom_name = None):  # optional : name in python
    call_params = [param(selfType, "self")] + parameters
    self.add_function_as_method("reference_as_pointer",
                                retval(ptr(returnType), reference_existing_object=True),
                                call_params,
                                template_parameters = [selfType, returnType, "&" + method],
                                foreign_cpp_namespace = "Spheral",
                                custom_name = custom_name)
    return

#-------------------------------------------------------------------------------
# const_ref_return_value
#-------------------------------------------------------------------------------
def const_ref_return_value(self,                 # The class we're adding to
                           selfType,             # The name of the class we're adding to
                           method,               # The accessor method returning a reference
                           returnType,           # The return type (without the &)
                           parameters = [],      # optional : list or tuple of paramters for the method.
                           custom_name = None):  # optional : name in python
    call_params = [param(selfType, "self")] + parameters
    self.add_function_as_method("const_reference_as_pointer",
                                retval(ptr(returnType), reference_existing_object=True),
                                call_params,
                                template_parameters = [selfType, returnType, "&" + method],
                                foreign_cpp_namespace = "Spheral",
                                custom_name = custom_name)
    return
