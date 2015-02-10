//---------------------------------Spheral++----------------------------------//
// PyFileIO -- A python friendly version of the FileIO interface, for use
// creating python FileIO objects.
//
// This class overrides the FileIO::read methods, since BPL has problems 
// handing back references to int, bool, and double.
//
// Created by JMO, Tue Dec 27 22:12:00 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_PyFileIO_hh__
#define __Spheral_PyFileIO_hh__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace FileIOSpace {

class PyFileIO: public FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  PyFileIO();
  PyFileIO(const std::string filename, AccessType access);

  // Destructor.
  virtual ~PyFileIO();

  // Descendent python objects should provide these write methods.
  virtual void write_unsigned_int(const unsigned value, const std::string pathName) = 0;
  virtual void write_int(const int value, const std::string pathName) = 0;
  virtual void write_bool(const bool value, const std::string pathName) = 0;
  virtual void write_double(const double value, const std::string pathName) = 0;
  virtual void write_string(const std::string value, const std::string pathName) = 0;

  virtual void write_Vector1d(const Dim<1>::Vector& value, const std::string pathName) = 0;
  virtual void write_Tensor1d(const Dim<1>::Tensor& value, const std::string pathName) = 0;
  virtual void write_SymTensor1d(const Dim<1>::SymTensor& value, const std::string pathName) = 0;
  virtual void write_ThirdRankTensor1d(const Dim<1>::ThirdRankTensor& value, const std::string pathName) = 0;

  virtual void write_Vector2d(const Dim<2>::Vector& value, const std::string pathName) = 0;
  virtual void write_Tensor2d(const Dim<2>::Tensor& value, const std::string pathName) = 0;
  virtual void write_SymTensor2d(const Dim<2>::SymTensor& value, const std::string pathName) = 0;
  virtual void write_ThirdRankTensor2d(const Dim<2>::ThirdRankTensor& value, const std::string pathName) = 0;

  virtual void write_Vector3d(const Dim<3>::Vector& value, const std::string pathName) = 0;
  virtual void write_Tensor3d(const Dim<3>::Tensor& value, const std::string pathName) = 0;
  virtual void write_SymTensor3d(const Dim<3>::SymTensor& value, const std::string pathName) = 0;
  virtual void write_ThirdRankTensor3d(const Dim<3>::ThirdRankTensor& value, const std::string pathName) = 0;

  virtual void write_vector_of_int(const std::vector<int>& value, const std::string pathName) = 0;
  virtual void write_vector_of_double(const std::vector<double>& value, const std::string pathName) = 0;
  virtual void write_vector_of_string(const std::vector<std::string>& value, const std::string pathName) = 0;

  virtual void write_vector_of_Vector1d(const std::vector<Dim<1>::Vector>& value, const std::string pathName) = 0;
  virtual void write_vector_of_Tensor1d(const std::vector<Dim<1>::Tensor>& value, const std::string pathName) = 0;
  virtual void write_vector_of_SymTensor1d(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName) = 0;
  virtual void write_vector_of_ThirdRankTensor1d(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) = 0;

  virtual void write_vector_of_Vector2d(const std::vector<Dim<2>::Vector>& value, const std::string pathName) = 0;
  virtual void write_vector_of_Tensor2d(const std::vector<Dim<2>::Tensor>& value, const std::string pathName) = 0;
  virtual void write_vector_of_SymTensor2d(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName) = 0;
  virtual void write_vector_of_ThirdRankTensor2d(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) = 0;

  virtual void write_vector_of_Vector3d(const std::vector<Dim<3>::Vector>& value, const std::string pathName) = 0;
  virtual void write_vector_of_Tensor3d(const std::vector<Dim<3>::Tensor>& value, const std::string pathName) = 0;
  virtual void write_vector_of_SymTensor3d(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName) = 0;
  virtual void write_vector_of_ThirdRankTensor3d(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) = 0;

  virtual void write_ScalarField1d(const FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) = 0;
  virtual void write_VectorField1d(const FieldSpace::Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) = 0;
  virtual void write_TensorField1d(const FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) = 0;
  virtual void write_SymTensorField1d(const FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write_ThirdRankTensorField1d(const FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write_IntField1d(const FieldSpace::Field<Dim<1>, int>& field, const std::string pathName) = 0;

  virtual void write_ScalarField2d(const FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) = 0;
  virtual void write_VectorField2d(const FieldSpace::Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) = 0;
  virtual void write_TensorField2d(const FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) = 0;
  virtual void write_SymTensorField2d(const FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write_ThirdRankTensorField2d(const FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write_IntField2d(const FieldSpace::Field<Dim<2>, int>& field, const std::string pathName) = 0;

  virtual void write_ScalarField3d(const FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) = 0;
  virtual void write_VectorField3d(const FieldSpace::Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) = 0;
  virtual void write_TensorField3d(const FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) = 0;
  virtual void write_SymTensorField3d(const FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write_ThirdRankTensorField3d(const FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write_IntField3d(const FieldSpace::Field<Dim<3>, int>& field, const std::string pathName) = 0;

  // Descendent python objects should provide these read methods.
  virtual void read_Vector1d(Dim<1>::Vector& value, const std::string pathName) const = 0;
  virtual void read_Tensor1d(Dim<1>::Tensor& value, const std::string pathName) const = 0;
  virtual void read_SymTensor1d(Dim<1>::SymTensor& value, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensor1d(Dim<1>::ThirdRankTensor& value, const std::string pathName) const = 0;

  virtual void read_Vector2d(Dim<2>::Vector& value, const std::string pathName) const = 0;
  virtual void read_Tensor2d(Dim<2>::Tensor& value, const std::string pathName) const = 0;
  virtual void read_SymTensor2d(Dim<2>::SymTensor& value, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensor2d(Dim<2>::ThirdRankTensor& value, const std::string pathName) const = 0;

  virtual void read_Vector3d(Dim<3>::Vector& value, const std::string pathName) const = 0;
  virtual void read_Tensor3d(Dim<3>::Tensor& value, const std::string pathName) const = 0;
  virtual void read_SymTensor3d(Dim<3>::SymTensor& value, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensor3d(Dim<3>::ThirdRankTensor& value, const std::string pathName) const = 0;

  virtual void read_vector_of_int(std::vector<int>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_double(std::vector<double>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_string(std::vector<std::string>& value, const std::string pathName) const = 0;

  virtual void read_vector_of_Vector1d(std::vector<Dim<1>::Vector>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_Tensor1d(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_SymTensor1d(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_ThirdRankTensor1d(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const = 0;

  virtual void read_vector_of_Vector2d(std::vector<Dim<2>::Vector>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_Tensor2d(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_SymTensor2d(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_ThirdRankTensor2d(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const = 0;

  virtual void read_vector_of_Vector3d(std::vector<Dim<3>::Vector>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_Tensor3d(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_SymTensor3d(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const = 0;
  virtual void read_vector_of_ThirdRankTensor3d(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const = 0;

  virtual void read_ScalarField1d(FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) const = 0;
  virtual void read_VectorField1d(FieldSpace::Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) const = 0;
  virtual void read_TensorField1d(FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) const = 0;
  virtual void read_SymTensorField1d(FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensorField1d(FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) const = 0;
  virtual void read_IntField1d(FieldSpace::Field<Dim<1>, int>& field, const std::string pathName) const = 0;

  virtual void read_ScalarField2d(FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) const = 0;
  virtual void read_VectorField2d(FieldSpace::Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) const = 0;
  virtual void read_TensorField2d(FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) const = 0;
  virtual void read_SymTensorField2d(FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensorField2d(FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) const = 0;
  virtual void read_IntField2d(FieldSpace::Field<Dim<2>, int>& field, const std::string pathName) const = 0;

  virtual void read_ScalarField3d(FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) const = 0;
  virtual void read_VectorField3d(FieldSpace::Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) const = 0;
  virtual void read_TensorField3d(FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) const = 0;
  virtual void read_SymTensorField3d(FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensorField3d(FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) const = 0;
  virtual void read_IntField3d(FieldSpace::Field<Dim<3>, int>& field, const std::string pathName) const = 0;

  // Override the base FileIO read methods to use the above pure virtual methods.

  // Write methods.
  virtual void write(const unsigned value, const std::string pathName) { write_unsigned_int(value, pathName); }
  virtual void write(const int value, const std::string pathName) { write_int(value, pathName); }
  virtual void write(const bool value, const std::string pathName) { write_bool(value, pathName); }
  virtual void write(const double value, const std::string pathName) { write_double(value, pathName); }
  virtual void write(const std::string value, const std::string pathName) { write_string(value, pathName); }

  virtual void write(const Dim<1>::Vector& value, const std::string pathName) { write_Vector1d(value, pathName); }
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName) { write_Tensor1d(value, pathName); }
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName) { write_SymTensor1d(value, pathName); }
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName) { write_ThirdRankTensor1d(value, pathName); }

  virtual void write(const Dim<2>::Vector& value, const std::string pathName) { write_Vector2d(value, pathName); }
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName) { write_Tensor2d(value, pathName); }
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName) { write_SymTensor2d(value, pathName); }
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName) { write_ThirdRankTensor2d(value, pathName); }

  virtual void write(const Dim<3>::Vector& value, const std::string pathName) { write_Vector3d(value, pathName); }
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName) { write_Tensor3d(value, pathName); }
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName) { write_SymTensor3d(value, pathName); }
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName) { write_ThirdRankTensor3d(value, pathName); }

  virtual void write(const std::vector<int>& value, const std::string pathName) { write_vector_of_int(value, pathName); }
  virtual void write(const std::vector<double>& value, const std::string pathName) { write_vector_of_double(value, pathName); }
  virtual void write(const std::vector<std::string>& value, const std::string pathName) { write_vector_of_string(value, pathName); }

  virtual void write(const std::vector<Dim<1>::Vector>& value, const std::string pathName) { write_vector_of_Vector1d(value, pathName); }
  virtual void write(const std::vector<Dim<1>::Tensor>& value, const std::string pathName) { write_vector_of_Tensor1d(value, pathName); }
  virtual void write(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName) { write_vector_of_SymTensor1d(value, pathName); }
  virtual void write(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) { write_vector_of_ThirdRankTensor1d(value, pathName); }

  virtual void write(const std::vector<Dim<2>::Vector>& value, const std::string pathName) { write_vector_of_Vector2d(value, pathName); }
  virtual void write(const std::vector<Dim<2>::Tensor>& value, const std::string pathName) { write_vector_of_Tensor2d(value, pathName); }
  virtual void write(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName) { write_vector_of_SymTensor2d(value, pathName); }
  virtual void write(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) { write_vector_of_ThirdRankTensor2d(value, pathName); }

  virtual void write(const std::vector<Dim<3>::Vector>& value, const std::string pathName) { write_vector_of_Vector3d(value, pathName); }
  virtual void write(const std::vector<Dim<3>::Tensor>& value, const std::string pathName) { write_vector_of_Tensor3d(value, pathName); }
  virtual void write(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName) { write_vector_of_SymTensor3d(value, pathName); }
  virtual void write(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) { write_vector_of_ThirdRankTensor3d(value, pathName); }

  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) { write_ScalarField1d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) { write_VectorField1d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) { write_TensorField1d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) { write_SymTensorField1d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) { write_ThirdRankTensorField1d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<1>, int>& value, const std::string pathName) { write_IntField1d(value, pathName); }

  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) { write_ScalarField2d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) { write_VectorField2d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) { write_TensorField2d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) { write_SymTensorField2d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) { write_ThirdRankTensorField2d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<2>, int>& value, const std::string pathName) { write_IntField2d(value, pathName); }

  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) { write_ScalarField3d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) { write_VectorField3d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) { write_TensorField3d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) { write_SymTensorField3d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) { write_ThirdRankTensorField3d(value, pathName); }
  virtual void write(const FieldSpace::Field<Dim<3>, int>& value, const std::string pathName) { write_IntField3d(value, pathName); }

  // Read methods.
  virtual void read(unsigned& value, const std::string pathName) const { value = read_unsigned_int(pathName); }
  virtual void read(int& value, const std::string pathName) const { value = read_int(pathName); }
  virtual void read(bool& value, const std::string pathName) const { value = read_bool(pathName); }
  virtual void read(double& value, const std::string pathName) const {value = read_double(pathName); }
  virtual void read(std::string& value, const std::string pathName) const { value = read_string(pathName); }

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const { read_Vector1d(value, pathName); }
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const { read_Tensor1d(value, pathName); }
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const { read_SymTensor1d(value, pathName); }
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const { read_ThirdRankTensor1d(value, pathName); }

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const { read_Vector2d(value, pathName); }
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const { read_Tensor2d(value, pathName); }
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const { read_SymTensor2d(value, pathName); }
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const { read_ThirdRankTensor2d(value, pathName); }

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const { read_Vector3d(value, pathName); }
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const { read_Tensor3d(value, pathName); }
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const { read_SymTensor3d(value, pathName); }
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const { read_ThirdRankTensor3d(value, pathName); }

  virtual void read(std::vector<int>& value, const std::string pathName) const { read_vector_of_int(value, pathName); }
  virtual void read(std::vector<double>& value, const std::string pathName) const { read_vector_of_double(value, pathName); }
  virtual void read(std::vector<std::string>& value, const std::string pathName) const { read_vector_of_string(value, pathName); }

  virtual void read(std::vector<Dim<1>::Vector>& value, const std::string pathName) const { read_vector_of_Vector1d(value, pathName); }
  virtual void read(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const { read_vector_of_Tensor1d(value, pathName); }
  virtual void read(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const { read_vector_of_SymTensor1d(value, pathName); }
  virtual void read(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const { read_vector_of_ThirdRankTensor1d(value, pathName); }

  virtual void read(std::vector<Dim<2>::Vector>& value, const std::string pathName) const { read_vector_of_Vector2d(value, pathName); }
  virtual void read(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const { read_vector_of_Tensor2d(value, pathName); }
  virtual void read(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const { read_vector_of_SymTensor2d(value, pathName); }
  virtual void read(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const { read_vector_of_ThirdRankTensor2d(value, pathName); }

  virtual void read(std::vector<Dim<3>::Vector>& value, const std::string pathName) const { read_vector_of_Vector3d(value, pathName); }
  virtual void read(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const { read_vector_of_Tensor3d(value, pathName); }
  virtual void read(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const { read_vector_of_SymTensor3d(value, pathName); }
  virtual void read(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const { read_vector_of_ThirdRankTensor3d(value, pathName); }

  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const { read_ScalarField1d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const { read_VectorField1d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const { read_TensorField1d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const { read_SymTensorField1d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const { read_ThirdRankTensorField1d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<1>, int>& value, const std::string pathName) const { read_IntField1d(value, pathName); }

  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const { read_ScalarField2d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const { read_VectorField2d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const { read_TensorField2d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const { read_SymTensorField2d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const { read_ThirdRankTensorField2d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<2>, int>& value, const std::string pathName) const { read_IntField2d(value, pathName); }

  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const { read_ScalarField3d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const { read_VectorField3d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const { read_TensorField3d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const { read_SymTensorField3d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const { read_ThirdRankTensorField3d(value, pathName); }
  virtual void read(FieldSpace::Field<Dim<3>, int>& value, const std::string pathName) const { read_IntField3d(value, pathName); }

private:
  //--------------------------- Private Interface ---------------------------//
  // Don't allow assignment.
  PyFileIO& operator=(const FileIO& rhs);
};

}
}

#endif
