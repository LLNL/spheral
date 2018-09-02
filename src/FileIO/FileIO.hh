//---------------------------------Spheral++----------------------------------//
// FileIO -- Provide the generic interface to file objects.
//
// Created by JMO, Tue Jul 11 22:19:37 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral__FileIO_hh__
#define __Spheral__FileIO_hh__

#ifndef CXXONLY
#include "Python.h"
#endif

#include "Geometry/Dimension.hh"

#include <vector>
#include <string>
#include <sstream>

namespace Spheral {

template<typename Dimension> class GeomPlane;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

// Define the standard file access types.
enum class AccessType {
  Undefined = -1,
  Create = 0,
  Read = 1,
  Write = 2,
  ReadWrite = 3
};

class FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  FileIO();
  FileIO(const std::string filename, AccessType access);

  // Destructor.
  virtual ~FileIO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const std::string fileName, AccessType access) = 0;
  virtual void close() = 0;

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const unsigned value, const std::string pathName) = 0;
  virtual void write(const int value, const std::string pathName) = 0;
  virtual void write(const bool value, const std::string pathName) = 0;
  virtual void write(const double value, const std::string pathName) = 0;
  virtual void write(const std::string value, const std::string pathName) = 0;

  virtual void write(const Dim<1>::Vector& value, const std::string pathName) = 0;
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName) = 0;
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName) = 0;
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName) = 0;

  virtual void write(const Dim<2>::Vector& value, const std::string pathName) = 0;
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName) = 0;
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName) = 0;
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName) = 0;

  virtual void write(const Dim<3>::Vector& value, const std::string pathName) = 0;
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName) = 0;
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName) = 0;
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName) = 0;

  virtual void read(unsigned& value, const std::string pathName) const = 0;
  virtual void read(int& value, const std::string pathName) const = 0;
  virtual void read(bool& value, const std::string pathName) const = 0;
  virtual void read(double& value, const std::string pathName) const = 0;
  virtual void read(std::string& value, const std::string pathName) const = 0;

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const = 0;
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const = 0;
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const = 0;
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const = 0;

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const = 0;
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const = 0;
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const = 0;
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const = 0;

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const = 0;
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const = 0;
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const = 0;
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const = 0;

  // We also require that FileIO objects write vectors of the primitive types.
  virtual void write(const std::vector<int>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<double>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<std::string>& value, const std::string pathName) = 0;

  virtual void write(const std::vector<Dim<1>::Vector>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<1>::Tensor>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) = 0;

  virtual void write(const std::vector<Dim<2>::Vector>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<2>::Tensor>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) = 0;

  virtual void write(const std::vector<Dim<3>::Vector>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<3>::Tensor>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName) = 0;
  virtual void write(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) = 0;

  virtual void read(std::vector<int>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<double>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<std::string>& value, const std::string pathName) const = 0;

  virtual void read(std::vector<Dim<1>::Vector>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const = 0;

  virtual void read(std::vector<Dim<2>::Vector>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const = 0;

  virtual void read(std::vector<Dim<3>::Vector>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const = 0;
  virtual void read(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const = 0;

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
#ifdef SPHERAL1D
  virtual void write(const Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<1>, int>& field, const std::string pathName) = 0;

  virtual void read(Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<1>, int>& field, const std::string pathName) const = 0;
#endif

#ifdef SPHERAL2D
  virtual void write(const Field<Dim<2>, Dim<1>::Scalar>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<2>, int>& field, const std::string pathName) = 0;

  virtual void read(Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<2>, int>& field, const std::string pathName) const = 0;
#endif

#ifdef SPHERAL3D
  virtual void write(const Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) = 0;
  virtual void write(const Field<Dim<3>, int>& field, const std::string pathName) = 0;

  virtual void read(Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) const = 0;
  virtual void read(Field<Dim<3>, int>& field, const std::string pathName) const = 0;
#endif
  //******************************************************************************

  // These methods are useful for the primitive types that are problematic
  // to return by reference from python.
  virtual void write_unsigned_int(const unsigned value, const std::string pathName) { this->write(value, pathName); }
  virtual void write_int(const int value, const std::string pathName)               { this->write(value, pathName); }
  virtual void write_bool(const bool value, const std::string pathName)             { this->write(value, pathName); }
  virtual void write_double(const double value, const std::string pathName)         { this->write(value, pathName); }
  virtual void write_string(const std::string value, const std::string pathName)    { this->write(value, pathName); }

  virtual unsigned read_unsigned_int(const std::string pathName) const { unsigned result;    this->read(result, pathName); return result; }
  virtual int read_int(const std::string pathName) const               { int result;         this->read(result, pathName); return result; }
  virtual bool read_bool(const std::string pathName) const             { bool result;        this->read(result, pathName); return result; }
  virtual double read_double(const std::string pathName) const         { double result;      this->read(result, pathName); return result; }
  virtual std::string read_string(const std::string pathName) const    { std::string result; this->read(result, pathName); return result; }

  // Read/write planes using the primitive methods.
  void write(const GeomPlane<Dim<1> >& value, const std::string pathName);
  void write(const GeomPlane<Dim<2> >& value, const std::string pathName);
  void write(const GeomPlane<Dim<3> >& value, const std::string pathName);

  void read(GeomPlane<Dim<1> >& value, const std::string pathName) const;
  void read(GeomPlane<Dim<2> >& value, const std::string pathName) const;
  void read(GeomPlane<Dim<3> >& value, const std::string pathName) const;

  // Provide char* read/write methods, which will simply call the string
  // methods provided by descendents.
  void write(const char* value, const std::string pathName);
  void read(char* value, const std::string pathName) const;

  // Provide float read/write methods, which will simply call the Scalar
  // methods provided by descendents.
  void write(const float& value, const std::string pathName);
  void read(float& value, const std::string pathName) const;

  void write(const std::vector<float>& value, const std::string pathName);
  void read(std::vector<float>& value, const std::string pathName) const;

  // Write/read FieldLists.
  template<typename Dimension, typename DataType>
  void write(const FieldList<Dimension, DataType>& fieldList,
             const std::string pathName);
  template<typename Dimension, typename DataType>
  void read(FieldList<Dimension, DataType>& fieldList,
            const std::string pathName) const;

  // Write/read Fields of vectors.
  template<typename Dimension, typename DataType>
  void write(const Field<Dimension, std::vector<DataType> >& field,
             const std::string pathName);
  template<typename Dimension, typename DataType>
  void read(Field<Dimension, std::vector<DataType> >& field,
            const std::string pathName) const;

  // Write/read a vector<DataType> if DataType is a primitive we already know about.
  template<typename DataType>
  void write(const std::vector<DataType>& x, const std::string pathName);
  template<typename DataType>
  void read(std::vector<DataType>& x, const std::string pathName) const;

  // A helper function to split a string up into substrings delimited by '/'.
  std::vector<std::string> splitPathComponents(const std::string path) const;

  // Return the group (directory) component of a path.
  std::string groupName(const std::string pathName) const;

  // Return the variable component of a path.
  std::string variableName(const std::string pathName) const;

  // Allow const access to some of the member data.
  const std::string& fileName() const;
  AccessType access() const;
  bool fileOpen() const;

#ifndef CXXONLY
  // These methods are particular to Python file objects.
  void writeObject(PyObject* thing, PyObject* path);
  PyObject* readObject(PyObject* path) const;
#endif

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendent class are allowed to directly diddle this common data.
  std::string mFileName;
  AccessType mAccess;
  bool mFileOpen;

#ifndef CXXONLY
  PyObject* mPickleMod;
  PyObject* mPickleDumps;
  PyObject* mPickleLoads;
#endif

private:
  //--------------------------- Private Interface ---------------------------//
  // Don't allow assignment.
  FileIO& operator=(const FileIO& rhs);
};

}

#include "FileIOInline.hh"

#else

// Forward declaration.
namespace Spheral{
  class FileIO;
}

#endif
