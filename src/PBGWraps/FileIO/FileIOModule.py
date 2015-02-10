from pybindgen import *

from PBGutils import *

#-------------------------------------------------------------------------------
# The class to handle wrapping this module.
#-------------------------------------------------------------------------------
class FileIO:

    FileIOTypes = [
        "Vector1d", "Tensor1d", "SymTensor1d", "ThirdRankTensor1d",
        "Vector2d", "Tensor2d", "SymTensor2d", "ThirdRankTensor2d",
        "Vector3d", "Tensor3d", "SymTensor3d", "ThirdRankTensor3d",
        "vector_of_int", "vector_of_double", "vector_of_string",
        "vector_of_Vector1d", "vector_of_Tensor1d", "vector_of_SymTensor1d", "vector_of_ThirdRankTensor1d",
        "vector_of_Vector2d", "vector_of_Tensor2d", "vector_of_SymTensor2d", "vector_of_ThirdRankTensor2d",
        "vector_of_Vector3d", "vector_of_Tensor3d", "vector_of_SymTensor3d", "vector_of_ThirdRankTensor3d",
        "Spheral::FieldSpace::ScalarField1d",
        "Spheral::FieldSpace::VectorField1d",
        "Spheral::FieldSpace::TensorField1d",
        "Spheral::FieldSpace::SymTensorField1d",
        "Spheral::FieldSpace::ThirdRankTensorField1d",
        "Spheral::FieldSpace::IntField1d",
        "Spheral::FieldSpace::ScalarField2d",
        "Spheral::FieldSpace::VectorField2d",
        "Spheral::FieldSpace::TensorField2d",
        "Spheral::FieldSpace::SymTensorField2d",
        "Spheral::FieldSpace::ThirdRankTensorField2d",
        "Spheral::FieldSpace::IntField2d",
        "Spheral::FieldSpace::ScalarField3d",
        "Spheral::FieldSpace::VectorField3d",
        "Spheral::FieldSpace::TensorField3d",
        "Spheral::FieldSpace::SymTensorField3d",
        "Spheral::FieldSpace::ThirdRankTensorField3d",
        "Spheral::FieldSpace::IntField3d",
        "double", "std::string", "int", "bool", "unsigned int",
        ]

    #---------------------------------------------------------------------------
    # Add the types to the given module.
    #---------------------------------------------------------------------------
    def __init__(self, mod, srcdir, topsrcdir):

        # Includes.
        mod.add_include('"%s/FileIO/FileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/FlatFileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/SiloFileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/PyFileIO.hh"' % topsrcdir)
        mod.add_include('"%s/FileIO/vectorstringUtilities.hh"' % topsrcdir)

        # Namespace.
        Spheral = mod.add_cpp_namespace("Spheral")
        space = Spheral.add_cpp_namespace("FileIOSpace")

        # Expose types.
        self.FileIO = addObject(space, "FileIO", allow_subclassing=True)
        self.FlatFileIO = addObject(space, "FlatFileIO", parent=self.FileIO, allow_subclassing=True)
        self.SiloFileIO = addObject(space, "SiloFileIO", parent=self.FileIO, allow_subclassing=True)
        self.PyFileIO = addObject(space, "PyFileIO", parent=self.FileIO, allow_subclassing=True)
        
        self.AccessType = space.add_enum("AccessType", 
                                         ["Undefined",
                                          "Create",
                                          "Read",
                                          "Write",
                                          "ReadWrite"])

        self.FlatFileFormat = space.add_enum("FlatFileFormat", ["ascii", "binary"])

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
    # The new sub modules (namespaces) introduced.
    #---------------------------------------------------------------------------
    def newSubModules(self):
        return ["FileIOSpace"]

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
                           param("FlatFileFormat", "format", default_value="Spheral::FileIOSpace::ascii")])


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

        # Add the read/write methods.
        for val in self.FileIOTypes:
            self._addPyFileIOReadWriteMethods(x, val)

        # Add our overrides for the base methods.
        for val in self.FileIOTypes:
            self._addFileIOReadWriteMethods(x, val, False)

        return

    #---------------------------------------------------------------------------
    # Helper for adding the virtual read/write methods for the given type
    # to a FileIO object.
    #---------------------------------------------------------------------------
    def _addFileIOReadWriteMethods(self, fio_obj, val, pureVirtual=True):
        if val in ["unsigned int", "int", "bool", "double", "std::string"]:
            fio_obj.add_method("write", None, [param(val, "value"),
                                               param("std::string", "pathName")],
                               is_virtual = True,
                               is_pure_virtual = pureVirtual)
            fio_obj.add_method("read_%s" % val.replace("std::", "").replace(" ", "_"),
                               val,
                               [param("std::string", "pathName")],
                               is_const = True,
                               is_virtual = True)
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
    # Helper for adding the virtual read/write methods for the given type
    # to a PyFileIO object.
    #---------------------------------------------------------------------------
    def _addPyFileIOReadWriteMethods(self, pyfio_obj, val):
        stripval = val.split("::")[-1]
        if val in ["unsigned int", "int", "bool", "double", "std::string"]:
            pyfio_obj.add_method("write_%s" % stripval.replace(" ", "_"),
                                 None,
                                 [param(val, "value"), 
                                  param("const std::string", "pathName")],
                                 is_virtual = True,
                                 is_pure_virtual = True)
            # pyfio_obj.add_method("read_%s" % stripval.replace(" ", "_"),
            #                      val,
            #                      [param("const std::string", "pathName")],
            #                      is_const = True,
            #                      is_virtual = True,
            #                      is_pure_virtual = True)
        else:
            pyfio_obj.add_method("write_%s" % stripval,
                                 None,
                                 [constrefparam(val, "value"),
                                  param("const std::string", "pathName")],
                                 is_virtual = True,
                                 is_pure_virtual = True)
            pyfio_obj.add_method("read_%s" % stripval,
                                 None, 
                                 [refparam(val, "value"),
                                  param("const std::string", "pathName")],
                                 is_const = True,
                                 is_virtual = True,
                                 is_pure_virtual = True)
        return
