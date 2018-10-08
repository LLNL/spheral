#-------------------------------------------------------------------------------
# PyFileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIO import *
from FileIOAbstractMethods import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11module("SpheralFileIO")
class PyFileIO(FileIO):
    "PyFileIO -- A python friendly version of the FileIO interface, for use creating python FileIO objects."

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Construct a PyFileIO object with a file name and access"

    #...........................................................................
    # Abstract interface
    @PYB11pure_virtual
    def write_Vector1d(self,
                       value = "const Dim<1>::Vector&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_Tensor1d(self,
                       value = "const Dim<1>::Tensor&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_SymTensor1d(self,
                          value = "const Dim<1>::SymTensor&",
                          pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_ThirdRankTensor1d(self,
                                value = "const Dim<1>::ThirdRankTensor&",
                                pathName = "const std::string"):
        return "void"

    @PYB11pure_virtual
    def write_Vector2d(self,
                       value = "const Dim<2>::Vector&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_Tensor2d(self,
                       value = "const Dim<2>::Tensor&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_SymTensor2d(self,
                          value = "const Dim<2>::SymTensor&",
                          pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_ThirdRankTensor2d(self,
                                value = "const Dim<2>::ThirdRankTensor&",
                                pathName = "const std::string"):
        return "void"

    @PYB11pure_virtual
    def write_Vector3d(self,
                       value = "const Dim<3>::Vector&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_Tensor3d(self,
                       value = "const Dim<3>::Tensor&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_SymTensor3d(self,
                          value = "const Dim<3>::SymTensor&",
                          pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_ThirdRankTensor3d(self,
                                value = "const Dim<3>::ThirdRankTensor&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_int(self,
                            value = "const std::vector<int>&",
                            pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_double(self,
                               value = "const std::vector<double>&",
                               pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_string(self,
                               value = "const std::vector<std::string>&",
                               pathName = "const std::string"):
        return "void"

    @PYB11pure_virtual
    def write_vector_of_Vector1d(self,
                                 value = "const std::vector<Dim<1>::Vector>&",
                                 pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_Tensor1d(self,
                                 value = "const std::vector<Dim<1>::Tensor>&",
                                 pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_SymTensor1d(self,
                                    value = "const std::vector<Dim<1>::SymTensor>&",
                                    pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_ThirdRankTensor1d(self,
                                          value = "const std::vector<Dim<1>::ThirdRankTensor>&",
                                          pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_Vector2d(self,
                                 value = "const std::vector<Dim<2>::Vector>&",
                                 pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_Tensor2d(self,
                                 value = "const std::vector<Dim<2>::Tensor>&",
                                 pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_SymTensor2d(self,
                                    value = "const std::vector<Dim<2>::SymTensor>&",
                                    pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_ThirdRankTensor2d(self,
                                          value = "const std::vector<Dim<2>::ThirdRankTensor>&",
                                          pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_Vector3d(self,
                                 value = "const std::vector<Dim<3>::Vector>&",
                                 pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_Tensor3d(self,
                                 value = "const std::vector<Dim<3>::Tensor>&",
                                 pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_SymTensor3d(self,
                                    value = "const std::vector<Dim<3>::SymTensor>&",
                                    pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    def write_vector_of_ThirdRankTensor3d(self,
                                          value = "const std::vector<Dim<3>::ThirdRankTensor>&",
                                          pathName = "const std::string"):
        return "void"

    @PYB11pure_virtual
    @PYB11const
    def read_Vector1d(self,
                      value = "Dim<1>::Vector&",
                      pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_Tensor1d(self,
                      value = "Dim<1>::Tensor&",
                      pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_SymTensor1d(self,
                         value = "Dim<1>::SymTensor&",
                         pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_ThirdRankTensor1d(self,
                               value = "Dim<1>::ThirdRankTensor&",
                               pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_Vector2d(self,
                      value = "Dim<2>::Vector&",
                      pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_Tensor2d(self,
                      value = "Dim<2>::Tensor&",
                      pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_SymTensor2d(self,
                         value = "Dim<2>::SymTensor&",
                         pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_ThirdRankTensor2d(self,
                               value = "Dim<2>::ThirdRankTensor&",
                               pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_Vector3d(self,
                      value = "Dim<3>::Vector&",
                      pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_Tensor3d(self,
                      value = "Dim<3>::Tensor&",
                      pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_SymTensor3d(self,
                         value = "Dim<3>::SymTensor&",
                         pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_ThirdRankTensor3d(self,
                               value = "Dim<3>::ThirdRankTensor&",
                               pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_int(self,
                           value = "std::vector<int>&",
                           pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_double(self,
                              value = "std::vector<double>&",
                              pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_string(self,
                              value = "std::vector<std::string>&",
                              pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_Vector1d(self,
                                value = "std::vector<Dim<1>::Vector>&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_Tensor1d(self,
                                value = "std::vector<Dim<1>::Tensor>&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_SymTensor1d(self,
                                   value = "std::vector<Dim<1>::SymTensor>&",
                                   pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_ThirdRankTensor1d(self,
                                         value = "std::vector<Dim<1>::ThirdRankTensor>&",
                                         pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_Vector2d(self,
                                value = "std::vector<Dim<2>::Vector>&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_Tensor2d(self,
                                value = "std::vector<Dim<2>::Tensor>&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_SymTensor2d(self,
                                   value = "std::vector<Dim<2>::SymTensor>&",
                                   pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_ThirdRankTensor2d(self,
                                         value = "std::vector<Dim<2>::ThirdRankTensor>&",
                                         pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_Vector3d(self,
                                value = "std::vector<Dim<3>::Vector>&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_Tensor3d(self,
                                value = "std::vector<Dim<3>::Tensor>&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_SymTensor3d(self,
                                   value = "std::vector<Dim<3>::SymTensor>&",
                                   pathName = "const std::string"):
        return "void"
  
    @PYB11pure_virtual
    @PYB11const
    def read_vector_of_ThirdRankTensor3d(self,
                                         value = "std::vector<Dim<3>::ThirdRankTensor>&",
                                         pathName = "const std::string"):
        return "void"

    for ndim in dims:
        exec('''
@PYB11pure_virtual
def write_ScalarField%(ndim)id(self,
                        field = "const Field<Dim<%(ndim)i>, Dim<%(ndim)i>::Scalar>&",
                        pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
def write_VectorField%(ndim)id(self,
                        field = "const Field<Dim<%(ndim)i>, Dim<%(ndim)i>::Vector>&",
                        pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
def write_TensorField%(ndim)id(self,
                        field = "const Field<Dim<%(ndim)i>, Dim<%(ndim)i>::Tensor>&",
                        pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
def write_SymTensorField%(ndim)id(self,
                           field = "const Field<Dim<%(ndim)i>, Dim<%(ndim)i>::SymTensor>&",
                           pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
def write_ThirdRankTensorField%(ndim)id(self,
                                 field = "const Field<Dim<%(ndim)i>, Dim<%(ndim)i>::ThirdRankTensor>&",
                                 pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
def write_IntField%(ndim)id(self,
                     field = "const Field<Dim<%(ndim)i>, int>&",
                     pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
@PYB11const
def read_ScalarField%(ndim)id(self,
                       field = "Field<Dim<%(ndim)i>, Dim<%(ndim)i>::Scalar>&",
                       pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
@PYB11const
def read_VectorField%(ndim)id(self,
                       field = "Field<Dim<%(ndim)i>, Dim<%(ndim)i>::Vector>&",
                       pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
@PYB11const
def read_TensorField%(ndim)id(self,
                       field = "Field<Dim<%(ndim)i>, Dim<%(ndim)i>::Tensor>&",
                       pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
@PYB11const
def read_SymTensorField%(ndim)id(self,
                          field = "Field<Dim<%(ndim)i>, Dim<%(ndim)i>::SymTensor>&",
                          pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
@PYB11const
def read_ThirdRankTensorField%(ndim)id(self,
                                field = "Field<Dim<%(ndim)i>, Dim<%(ndim)i>::ThirdRankTensor>&",
                                pathName = "const std::string"):
    return "void"

@PYB11pure_virtual
@PYB11const
def read_IntField%(ndim)id(self,
                    field = "Field<Dim<%(ndim)i>, int>&",
                    pathName = "const std::string"):
    return "void"
''' % {"ndim" : ndim})

#-------------------------------------------------------------------------------
# Override the required virtual interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, PyFileIO, virtual=True, pure_virtual=False)
