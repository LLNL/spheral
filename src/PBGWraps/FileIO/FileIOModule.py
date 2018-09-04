from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class FileIO:

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir, dims):

        self.dims = dims

        self.FileIOTypes = []
        self.FileIOTemplateTypes = []
        for ndim in (1, 2, 3):
            exec("""
self.FileIOTypes += [
        "Vector%(dim)s", "Tensor%(dim)s", "SymTensor%(dim)s", "ThirdRankTensor%(dim)s",
        "vector_of_Vector%(dim)s",
        "vector_of_Tensor%(dim)s",
        "vector_of_SymTensor%(dim)s",
        "vector_of_ThirdRankTensor%(dim)s",
        ]""" % {"dim" : "%id" % ndim,
                "Dim" : "Spheral::Dim<%i>" % ndim})
        for ndim in self.dims:
            exec("""
self.FileIOTypes += [
        "Spheral::ScalarField%(dim)s",
        "Spheral::VectorField%(dim)s",
        "Spheral::TensorField%(dim)s",
        "Spheral::SymTensorField%(dim)s",
        "Spheral::ThirdRankTensorField%(dim)s",
        "Spheral::IntField%(dim)s",
        ]
self.FileIOTemplateTypes += [
        ("Spheral::ScalarFieldList%(dim)s", ["%(Dim)s", "double"]),
        ("Spheral::VectorFieldList%(dim)s", ["%(Dim)s", "Vector%(dim)s"]),
        ("Spheral::TensorFieldList%(dim)s", ["%(Dim)s", "Tensor%(dim)s"]),
        ("Spheral::SymTensorFieldList%(dim)s", ["%(Dim)s", "SymTensor%(dim)s"]),
        ("Spheral::ThirdRankTensorFieldList%(dim)s", ["%(Dim)s", "ThirdRankTensor%(dim)s"]),
        ("Spheral::IntFieldList%(dim)s", ["%(Dim)s", "int"]),
        ("Spheral::VectorDoubleFieldList%(dim)s", ["%(Dim)s", "vector_of_double"]),
        ]""" % {"dim" : "%id" % ndim,
                "Dim" : "Spheral::Dim<%i>" % ndim})
        self.FileIOTypes += ["vector_of_int", "vector_of_double", "vector_of_string",
                             "double", "std::string", "int", "bool", "unsigned int"]

        # Includes.
        mod.add_include('"%s/FileIO/FileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/FlatFileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/SiloFileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/PyFileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/vectorstringUtilities.hh"' % topsrcdir)

        # Namespace.
        space = mod.add_cpp_namespace("Spheral")

        # Expose types.
        self.FileIO = addObject(space, "FileIO", allow_subclassing=True)
        self.FlatFileIO = addObject(space, "FlatFileIO", parent=self.FileIO, allow_subclassing=True)
        self.SiloFileIO = addObject(space, "SiloFileIO", parent=self.FileIO, allow_subclassing=True)
        self.PyFileIO = addObject(space, "PyFileIO", parent=self.FileIO, allow_subclassing=True)
        
        self.AccessType = space.add_enum("AccessType", [("Undefined", "Spheral::AccessType::Undefined"),
                                                        ("Create", "Spheral::AccessType::Create"),
                                                        ("Read", "Spheral::AccessType::Read"),
                                                        ("Write", "Spheral::AccessType::Write"),
                                                        ("ReadWrite", "Spheral::AccessType::ReadWrite")])

        self.FlatFileFormat = space.add_enum("FlatFileFormat", [("ascii", "Spheral::FlatFileFormat::ascii"),
                                                                ("binary", "Spheral::FlatFileFormat::binary")])

        return

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def generateBindings(self, mod):

        self.addFileIOMethods()
        self.addFlatFileIOMethods()
        self.addSiloFileIOMethods()
        self.addPyFileIOMethods()

        # Add the functions.
        Spheral = mod.add_cpp_namespace("Spheral")
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_int", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_ULL", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_double", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_string", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_Vector1d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_Vector2d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_Vector3d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_Tensor1d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_Tensor2d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_Tensor3d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_SymTensor1d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_SymTensor2d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_SymTensor3d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_ThirdRankTensor1d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_ThirdRankTensor2d", "val"), param("int", "precision", default_value="30")])
        Spheral.add_function("vector2string", "std::string", [constrefparam("vector_of_ThirdRankTensor3d", "val"), param("int", "precision", default_value="30")])

        Spheral.add_function("string2vector", "vector_of_int", [param("std::string", "val")], template_parameters=["int"], custom_name="string2vector_of_int")
        Spheral.add_function("string2vector", "vector_of_unsigned", [param("std::string", "val")], template_parameters=["unsigned"], custom_name="string2vector_of_unsigned")
        Spheral.add_function("string2vector", "vector_of_ULL", [param("std::string", "val")], template_parameters=["uint64_t"], custom_name="string2vector_of_ULL")
        Spheral.add_function("string2vector", "vector_of_double", [param("std::string", "val")], template_parameters=["double"], custom_name="string2vector_of_double")
        Spheral.add_function("string2vector", "vector_of_string", [param("std::string", "val")], template_parameters=["std::string"], custom_name="string2vector_of_string")
        Spheral.add_function("string2vector", "vector_of_Vector1d", [param("std::string", "val")], template_parameters=["Vector1d"], custom_name="string2vector_of_Vector1d")
        Spheral.add_function("string2vector", "vector_of_Vector2d", [param("std::string", "val")], template_parameters=["Vector2d"], custom_name="string2vector_of_Vector2d")
        Spheral.add_function("string2vector", "vector_of_Vector3d", [param("std::string", "val")], template_parameters=["Vector3d"], custom_name="string2vector_of_Vector3d")
        Spheral.add_function("string2vector", "vector_of_Tensor1d", [param("std::string", "val")], template_parameters=["Tensor1d"], custom_name="string2vector_of_Tensor1d")
        Spheral.add_function("string2vector", "vector_of_Tensor2d", [param("std::string", "val")], template_parameters=["Tensor2d"], custom_name="string2vector_of_Tensor2d")
        Spheral.add_function("string2vector", "vector_of_Tensor3d", [param("std::string", "val")], template_parameters=["Tensor3d"], custom_name="string2vector_of_Tensor3d")
        Spheral.add_function("string2vector", "vector_of_SymTensor1d", [param("std::string", "val")], template_parameters=["SymTensor1d"], custom_name="string2vector_of_SymTensor1d")
        Spheral.add_function("string2vector", "vector_of_SymTensor2d", [param("std::string", "val")], template_parameters=["SymTensor2d"], custom_name="string2vector_of_SymTensor2d")
        Spheral.add_function("string2vector", "vector_of_SymTensor3d", [param("std::string", "val")], template_parameters=["SymTensor3d"], custom_name="string2vector_of_SymTensor3d")
        Spheral.add_function("string2vector", "vector_of_ThirdRankTensor1d", [param("std::string", "val")], template_parameters=["ThirdRankTensor1d"], custom_name="string2vector_of_ThirdRankTensor1d")
        Spheral.add_function("string2vector", "vector_of_ThirdRankTensor2d", [param("std::string", "val")], template_parameters=["ThirdRankTensor2d"], custom_name="string2vector_of_ThirdRankTensor2d")
        Spheral.add_function("string2vector", "vector_of_ThirdRankTensor3d", [param("std::string", "val")], template_parameters=["ThirdRankTensor3d"], custom_name="string2vector_of_ThirdRankTensor3d")

        return

    #---------------------------------------------------------------------------
    # Add FileIO methods.
    #---------------------------------------------------------------------------
    def addFileIOMethods(self):

        x = self.FileIO

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("std::string", "filename"),
                           param("AccessType", "access")])

        # Methods.
        x.add_method("open", None, [param("std::string", "name"),
                                    param("AccessType", "access")],
                     is_pure_virtual=True)
        x.add_method("close", None, [], is_pure_virtual=True)

        # Add the templated read/write methods.
        for val, template_params in self.FileIOTemplateTypes:
            self._addFileIOReadWriteTemplateMethods(x, val, template_params)

        # Add the standard read/write methods for the supported types.
        for val in self.FileIOTypes:
            self._addFileIOReadWriteMethods(x, val)

        # Write and read objects.
        x.add_method("writeObject", None, 
                     [param(ptr("PyObject"), "thing", transfer_ownership=False),
                      param(ptr("PyObject"), "path", transfer_ownership=False)])
        x.add_method("readObject", retval(ptr("PyObject"), caller_owns_return=True),
                     [param(ptr("PyObject"), "path", transfer_ownership=False)],
                     is_const = True)

        return

    #---------------------------------------------------------------------------
    # Add FlatFileIO methods.
    #---------------------------------------------------------------------------
    def addFlatFileIOMethods(self):

        x = self.FlatFileIO

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("std::string", "filename"),
                           param("AccessType", "access"),
                           param("FlatFileFormat", "format", default_value="Spheral::FlatFileFormat::ascii")])


        # Methods.
        x.add_method("open", None, [param("std::string", "name"),
                                    param("AccessType", "access")],
                     is_virtual=True)
        x.add_method("close", None, [], is_virtual=True)
        x.add_method("findPathName", None, [constrefparam("std::string", "pathName")], is_const=True)
        x.add_method("beginningOfFile", None, [], is_const=True)

        # Attributes.
        x.add_instance_attribute("precision", "int", getter="precision", setter="setPrecision")
        x.add_instance_attribute("readyToWrite", "bool", getter="readyToWrite", is_const=True)
        x.add_instance_attribute("readyToRead", "bool", getter="readyToRead", is_const=True)

        # Add the templated read/write methods.
        for val, template_params in self.FileIOTemplateTypes:
            self._addFileIOReadWriteTemplateMethods(x, val, template_params)

        # Add the standard read/write methods for the supported types.
        for val in self.FileIOTypes:
            self._addFileIOReadWriteMethods(x, val, False)

        return


    #---------------------------------------------------------------------------
    # Add SiloFileIO methods.
    #---------------------------------------------------------------------------
    def addSiloFileIOMethods(self):

        x = self.SiloFileIO

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("std::string", "filename"),
                           param("AccessType", "access")])

        # Methods.
        x.add_method("open", None, [param("std::string", "name"),
                                    param("AccessType", "access")],
                     is_virtual=True)
        x.add_method("close", None, [], is_virtual=True)

        # Add the templated read/write methods.
        for val, template_params in self.FileIOTemplateTypes:
            self._addFileIOReadWriteTemplateMethods(x, val, template_params)

        # Add the standard read/write methods for the supported types.
        for val in self.FileIOTypes:
            self._addFileIOReadWriteMethods(x, val, False)

        return

    #---------------------------------------------------------------------------
    # Add PyFileIO methods.
    #---------------------------------------------------------------------------
    def addPyFileIOMethods(self):

        x = self.PyFileIO

        # Constructors.
        x.add_constructor([])
        x.add_constructor([param("std::string", "filename"),
                           param("AccessType", "access")])

        # Add the templated read/write methods.
        for val, template_params in self.FileIOTemplateTypes:
            self._addFileIOReadWriteTemplateMethods(x, val, template_params)

        # Add our overrides for the base methods.
        for val in self.FileIOTypes:
            self._addFileIOReadWriteMethods(x, val, False)

        # Add the read/write methods.
        for val in self.FileIOTypes:
            self._addPyFileIOReadWriteMethods(x, val)

        return

    #---------------------------------------------------------------------------
    # Helper for adding the virtual read/write methods for the given type
    # to a FileIO object.
    #---------------------------------------------------------------------------
    def _addFileIOReadWriteMethods(self, fio_obj, val, pureVirtual=True):
        if val in ["unsigned int", "int", "bool", "double", "std::string"]:
            valname = val.replace("std::", "").replace(" ", "_")
            fio_obj.add_method("write_%s" % valname,
                               None,
                               [param(val, "value"), param("std::string", "pathName")],
                               is_virtual = True,
                               is_pure_virtual = pureVirtual)
            fio_obj.add_method("read_%s" % valname,
                               val,
                               [param("std::string", "pathName")],
                               is_const = True,
                               is_virtual = True)
            fio_obj.add_method("write", None, [param(val, "value"),
                                               param("std::string", "pathName")],
                               is_virtual = True,
                               is_pure_virtual = pureVirtual)
        else:
            fio_obj.add_method("write", None, [constrefparam(val, "value"),
                                               param("std::string", "pathName")],
                               is_virtual = True,
                               is_pure_virtual = pureVirtual)
        fio_obj.add_method("read", None, [refparam(val, "value"),
                                          param("std::string", "pathName")],
                           is_const = True,
                           is_virtual = True,
                           is_pure_virtual = pureVirtual)
        return

    #---------------------------------------------------------------------------
    # Helper for adding the templated read/write methods for the given type
    # to a FileIO object.
    #---------------------------------------------------------------------------
    def _addFileIOReadWriteTemplateMethods(self, fio_obj, val, template_params):
        fio_obj.add_method("write", None, [constrefparam(val, "value"),
                                           param("std::string", "pathName")],
                           template_parameters = template_params,
                           custom_name = "write")
        fio_obj.add_method("read", None, [refparam(val, "value"),
                                          param("std::string", "pathName")],
                           template_parameters = template_params,
                           is_const = True,
                           custom_name = "read")
        return

    #---------------------------------------------------------------------------
    # Helper for adding the virtual read/write methods for the given type
    # to a PyFileIO object.
    #---------------------------------------------------------------------------
    def _addPyFileIOReadWriteMethods(self, pyfio_obj, val):
        stripval = val.split("::")[-1]
        if (not val in ["unsigned int", "int", "bool", "double", "std::string"]):
            pyfio_obj.add_method("write_%s" % stripval,
                                 None,
                                 [constrefparam(val, "value"),
                                  param("const std::string", "pathName")],
                                 is_pure_virtual = True)
            pyfio_obj.add_method("read_%s" % stripval,
                                 None, 
                                 [refparam(val, "value"),
                                  param("const std::string", "pathName")],
                                 is_const = True,
                                 is_pure_virtual = True)
        return
