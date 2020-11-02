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

#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include <vector>

namespace Spheral {

class PyFileIO: public FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  PyFileIO();
  PyFileIO(const std::string filename, AccessType access);

  // Destructor.
  virtual ~PyFileIO();

  // Check if the specified path is in the file.
  virtual bool pathExists(const std::string) const override { VERIFY2(false, "pathExists not overridden"); }

  // Descendent python objects should provide these write methods.
  virtual void write_unsigned_int(const unsigned, const std::string)              override { VERIFY2(false, "write_unsigned_int not overridden"); }
  virtual void write_int(const int, const std::string)                            override { VERIFY2(false, "write_int not overridden"); }
  virtual void write_bool(const bool, const std::string)                          override { VERIFY2(false, "write_bool not overridden"); }
  virtual void write_double(const double, const std::string)                      override { VERIFY2(false, "write_double not overridden"); }
  virtual void write_string(const std::string, const std::string)                 override { VERIFY2(false, "write_string not overridden"); }
  virtual void write_vector_int(const std::vector<int>&, const std::string)                { VERIFY2(false, "write_vector_int not overridden"); }
  virtual void write_vector_double(const std::vector<double>&, const std::string)          { VERIFY2(false, "write_vector_double not overridden"); }
  virtual void write_vector_string(const std::vector<std::string>&, const std::string)     { VERIFY2(false, "write_vector_string not overridden"); }

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

#ifdef SPHERAL1D
  virtual void write_ScalarField1d(const Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) = 0;
  virtual void write_VectorField1d(const Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) = 0;
  virtual void write_TensorField1d(const Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) = 0;
  virtual void write_SymTensorField1d(const Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write_ThirdRankTensorField1d(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write_IntField1d(const Field<Dim<1>, int>& field, const std::string pathName) = 0;
#endif

#ifdef SPHERAL2D
  virtual void write_ScalarField2d(const Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) = 0;
  virtual void write_VectorField2d(const Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) = 0;
  virtual void write_TensorField2d(const Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) = 0;
  virtual void write_SymTensorField2d(const Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write_ThirdRankTensorField2d(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write_IntField2d(const Field<Dim<2>, int>& field, const std::string pathName) = 0;
#endif

#ifdef SPHERAL3D
  virtual void write_ScalarField3d(const Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) = 0;
  virtual void write_VectorField3d(const Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) = 0;
  virtual void write_TensorField3d(const Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) = 0;
  virtual void write_SymTensorField3d(const Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write_ThirdRankTensorField3d(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write_IntField3d(const Field<Dim<3>, int>& field, const std::string pathName) = 0;
#endif

  // Descendent python objects should provide these read methods.
  virtual unsigned read_unsigned_int(const std::string)                          const override { VERIFY2(false, "read_unsigned_int not overridden"); }
  virtual int read_int(const std::string)                                        const override { VERIFY2(false, "read_int not overridden"); }
  virtual bool read_bool(const std::string)                                      const override { VERIFY2(false, "read_bool not overridden"); }
  virtual double read_double(const std::string)                                  const override { VERIFY2(false, "read_double not overridden"); }
  virtual std::string read_string(const std::string)                             const override { VERIFY2(false, "read_string not overridden"); }
  virtual void read_vector_int(std::vector<int>*, const std::string) const                { VERIFY2(false, "read_vector_int not overridden"); }
  virtual void read_vector_double(std::vector<double>*, const std::string) const          { VERIFY2(false, "read_vector_double not overridden"); }
  virtual void read_vector_string(std::vector<std::string>*, const std::string) const     { VERIFY2(false, "read_vector_string not overridden"); }

  virtual Dim<1>::Vector read_Vector1d(const std::string pathName) const = 0;
  virtual Dim<1>::Tensor read_Tensor1d(const std::string pathName) const = 0;
  virtual Dim<1>::SymTensor read_SymTensor1d(const std::string pathName) const = 0;
  virtual Dim<1>::ThirdRankTensor read_ThirdRankTensor1d(const std::string pathName) const = 0;

  virtual Dim<2>::Vector read_Vector2d(const std::string pathName) const = 0;
  virtual Dim<2>::Tensor read_Tensor2d(const std::string pathName) const = 0;
  virtual Dim<2>::SymTensor read_SymTensor2d(const std::string pathName) const = 0;
  virtual Dim<2>::ThirdRankTensor read_ThirdRankTensor2d(const std::string pathName) const = 0;

  virtual Dim<3>::Vector read_Vector3d(const std::string pathName) const = 0;
  virtual Dim<3>::Tensor read_Tensor3d(const std::string pathName) const = 0;
  virtual Dim<3>::SymTensor read_SymTensor3d(const std::string pathName) const = 0;
  virtual Dim<3>::ThirdRankTensor read_ThirdRankTensor3d(const std::string pathName) const = 0;

#ifdef SPHERAL1D
  virtual void read_ScalarField1d(Field<Dim<1>, Dim<1>::Scalar>* field, const std::string pathName) const = 0;
  virtual void read_VectorField1d(Field<Dim<1>, Dim<1>::Vector>* field, const std::string pathName) const = 0;
  virtual void read_TensorField1d(Field<Dim<1>, Dim<1>::Tensor>* field, const std::string pathName) const = 0;
  virtual void read_SymTensorField1d(Field<Dim<1>, Dim<1>::SymTensor>* field, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensorField1d(Field<Dim<1>, Dim<1>::ThirdRankTensor>* field, const std::string pathName) const = 0;
  virtual void read_IntField1d(Field<Dim<1>, int>* field, const std::string pathName) const = 0;
#endif

#ifdef SPHERAL2D
  virtual void read_ScalarField2d(Field<Dim<2>, Dim<2>::Scalar>* field, const std::string pathName) const = 0;
  virtual void read_VectorField2d(Field<Dim<2>, Dim<2>::Vector>* field, const std::string pathName) const = 0;
  virtual void read_TensorField2d(Field<Dim<2>, Dim<2>::Tensor>* field, const std::string pathName) const = 0;
  virtual void read_SymTensorField2d(Field<Dim<2>, Dim<2>::SymTensor>* field, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensorField2d(Field<Dim<2>, Dim<2>::ThirdRankTensor>* field, const std::string pathName) const = 0;
  virtual void read_IntField2d(Field<Dim<2>, int>* field, const std::string pathName) const = 0;
#endif

#ifdef SPHERAL3D
  virtual void read_ScalarField3d(Field<Dim<3>, Dim<3>::Scalar>* field, const std::string pathName) const = 0;
  virtual void read_VectorField3d(Field<Dim<3>, Dim<3>::Vector>* field, const std::string pathName) const = 0;
  virtual void read_TensorField3d(Field<Dim<3>, Dim<3>::Tensor>* field, const std::string pathName) const = 0;
  virtual void read_SymTensorField3d(Field<Dim<3>, Dim<3>::SymTensor>* field, const std::string pathName) const = 0;
  virtual void read_ThirdRankTensorField3d(Field<Dim<3>, Dim<3>::ThirdRankTensor>* field, const std::string pathName) const = 0;
  virtual void read_IntField3d(Field<Dim<3>, int>* field, const std::string pathName) const = 0;
#endif

  //***************************************************************************
  // Override the base FileIO read methods to use the above pure virtual methods.
  // Write methods.
  virtual void write(const unsigned& value, const std::string pathName) override                 { write_unsigned_int(value, pathName); }
  virtual void write(const int& value, const std::string pathName) override                      { write_int(value, pathName); }
  virtual void write(const bool& value, const std::string pathName) override                     { write_bool(value, pathName); }
  virtual void write(const double& value, const std::string pathName) override                   { write_double(value, pathName); }
  virtual void write(const std::string& value, const std::string pathName) override              { write_string(value, pathName); }
  virtual void write(const std::vector<int>& value, const std::string pathName) override         { write_vector_int(value, pathName); }
  virtual void write(const std::vector<double>& value, const std::string pathName) override      { write_vector_double(value, pathName); }
  virtual void write(const std::vector<std::string>& value, const std::string pathName) override { write_vector_string(value, pathName); }

  virtual void write(const Dim<1>::Vector& value, const std::string pathName) override { write_Vector1d(value, pathName); }
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName) override { write_Tensor1d(value, pathName); }
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName) override { write_SymTensor1d(value, pathName); }
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName) override { write_ThirdRankTensor1d(value, pathName); }

  virtual void write(const Dim<2>::Vector& value, const std::string pathName) override { write_Vector2d(value, pathName); }
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName) override { write_Tensor2d(value, pathName); }
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName) override { write_SymTensor2d(value, pathName); }
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName) override { write_ThirdRankTensor2d(value, pathName); }

  virtual void write(const Dim<3>::Vector& value, const std::string pathName) override { write_Vector3d(value, pathName); }
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName) override { write_Tensor3d(value, pathName); }
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName) override { write_SymTensor3d(value, pathName); }
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName) override { write_ThirdRankTensor3d(value, pathName); }

#ifdef SPHERAL1D
  virtual void write(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) override { write_ScalarField1d(value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) override { write_VectorField1d(value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) override { write_TensorField1d(value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) override { write_SymTensorField1d(value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { write_ThirdRankTensorField1d(value, pathName); }
  virtual void write(const Field<Dim<1>, int>& value, const std::string pathName) override { write_IntField1d(value, pathName); }
#endif

#ifdef SPHERAL2D
  virtual void write(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) override { write_ScalarField2d(value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) override { write_VectorField2d(value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) override { write_TensorField2d(value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) override { write_SymTensorField2d(value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { write_ThirdRankTensorField2d(value, pathName); }
  virtual void write(const Field<Dim<2>, int>& value, const std::string pathName) override { write_IntField2d(value, pathName); }
#endif

#ifdef SPHERAL3D
  virtual void write(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) override { write_ScalarField3d(value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) override { write_VectorField3d(value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) override { write_TensorField3d(value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) override { write_SymTensorField3d(value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { write_ThirdRankTensorField3d(value, pathName); }
  virtual void write(const Field<Dim<3>, int>& value, const std::string pathName) override { write_IntField3d(value, pathName); }
#endif

  // Read methods.
  virtual void read(unsigned& value, const std::string pathName) const override                 { value = read_unsigned_int(pathName); }
  virtual void read(int& value, const std::string pathName) const override                      { value = read_int(pathName); }
  virtual void read(bool& value, const std::string pathName) const override                     { value = read_bool(pathName); }
  virtual void read(double& value, const std::string pathName) const override                   { value = read_double(pathName); }
  virtual void read(std::string& value, const std::string pathName) const override              { value = read_string(pathName); }
  virtual void read(std::vector<int>& value, const std::string pathName) const override         { read_vector_int(&value, pathName); }
  virtual void read(std::vector<double>& value, const std::string pathName) const override      { read_vector_double(&value, pathName); }
  virtual void read(std::vector<std::string>& value, const std::string pathName) const override { read_vector_string(&value, pathName); }

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const override { value = read_Vector1d(pathName); }
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const override { value = read_Tensor1d(pathName); }
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const override { value = read_SymTensor1d(pathName); }
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const override { value = read_ThirdRankTensor1d(pathName); }

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const override { value = read_Vector2d(pathName); }
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const override { value = read_Tensor2d(pathName); }
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const override { value = read_SymTensor2d(pathName); }
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const override { value = read_ThirdRankTensor2d(pathName); }

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const override { value = read_Vector3d(pathName); }
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const override { value = read_Tensor3d(pathName); }
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const override { value = read_SymTensor3d(pathName); }
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const override { value = read_ThirdRankTensor3d(pathName); }

#ifdef SPHERAL1D
  virtual void read(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const override { read_ScalarField1d(&value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const override { read_VectorField1d(&value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const override { read_TensorField1d(&value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const override { read_SymTensorField1d(&value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { read_ThirdRankTensorField1d(&value, pathName); }
  virtual void read(Field<Dim<1>, int>& value, const std::string pathName) const override { read_IntField1d(&value, pathName); }
#endif

#ifdef SPHERAL2D
  virtual void read(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const override { read_ScalarField2d(&value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const override { read_VectorField2d(&value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const override { read_TensorField2d(&value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const override { read_SymTensorField2d(&value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { read_ThirdRankTensorField2d(&value, pathName); }
  virtual void read(Field<Dim<2>, int>& value, const std::string pathName) const override { read_IntField2d(&value, pathName); }
#endif

#ifdef SPHERAL3D
  virtual void read(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const override { read_ScalarField3d(&value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const override { read_VectorField3d(&value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const override { read_TensorField3d(&value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const override { read_SymTensorField3d(&value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { read_ThirdRankTensorField3d(&value, pathName); }
  virtual void read(Field<Dim<3>, int>& value, const std::string pathName) const override { read_IntField3d(&value, pathName); }
#endif
  //***************************************************************************

  //------------------------------------------------------------------------------
  // We have to forward the templated write/read methods to the base class due to
  // function hiding.
  // Write/read a vector<Value> if Value is a primitive we already know about.
  template<typename Value> void write(const std::vector<Value>& x, const std::string pathName) { FileIO::write(x, pathName); }
  template<typename Value> void  read(std::vector<Value>& x, const std::string pathName) const { FileIO::read(x, pathName); }
  //------------------------------------------------------------------------------

  // Forward hidden base class methods
  using FileIO::write;
  using FileIO::read;

private:
  //--------------------------- Private Interface ---------------------------//
  // Don't allow assignment.
  PyFileIO& operator=(const FileIO& rhs);
};

}

#endif
