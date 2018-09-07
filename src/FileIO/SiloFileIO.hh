//---------------------------------Spheral++----------------------------------//
// SiloFileIO -- Provide the interface to silo file objects.
//
// Created by JMO, Sat Feb  7 23:06:03 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_SiloFileIO__
#define __Spheral_SiloFileIO__

#include "FileIO.hh"

#include <vector>
#include <string>

extern "C" {
#include "silo.h"
}

namespace Spheral {

class SiloFileIO: public FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  SiloFileIO();
  SiloFileIO(const std::string fileName, AccessType access);

  // Destructor.
  virtual ~SiloFileIO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const std::string fileName, AccessType access);
  virtual void close();

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const unsigned value, const std::string pathName);
  virtual void write(const int value, const std::string pathName);
  virtual void write(const bool value, const std::string pathName);
  virtual void write(const double value, const std::string pathName);
  virtual void write(const std::string value, const std::string pathName);

  virtual void write(const Dim<1>::Vector& value, const std::string pathName);
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName);
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName);
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName);

  virtual void write(const Dim<2>::Vector& value, const std::string pathName);
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName);
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName);
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName);

  virtual void write(const Dim<3>::Vector& value, const std::string pathName);
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName);
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName);
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName);

  virtual void read(unsigned& value, const std::string pathName) const;
  virtual void read(int& value, const std::string pathName) const;
  virtual void read(bool& value, const std::string pathName) const;
  virtual void read(double& value, const std::string pathName) const;
  virtual void read(std::string& value, const std::string pathName) const;

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const;
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const;
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const;
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const;

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const;
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const;
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const;
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const;

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const;
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const;
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const;
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const;

  // We also require that FileIO objects write vectors of the primitive types.
  virtual void write(const std::vector<int>& value, const std::string pathName);
  virtual void write(const std::vector<double>& value, const std::string pathName);
  virtual void write(const std::vector<std::string>& value, const std::string pathName);

  virtual void write(const std::vector<Dim<1>::Vector>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<1>::Tensor>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName);

  virtual void write(const std::vector<Dim<2>::Vector>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<2>::Tensor>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName);

  virtual void write(const std::vector<Dim<3>::Vector>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<3>::Tensor>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName);
  virtual void write(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName);

  virtual void read(std::vector<int>& value, const std::string pathName) const;
  virtual void read(std::vector<double>& value, const std::string pathName) const;
  virtual void read(std::vector<std::string>& value, const std::string pathName) const;

  virtual void read(std::vector<Dim<1>::Vector>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const;

  virtual void read(std::vector<Dim<2>::Vector>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const;

  virtual void read(std::vector<Dim<3>::Vector>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const;
  virtual void read(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const;

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
#ifdef SPHERAL1D
  virtual void write(const Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName);
  virtual void write(const Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName);
  virtual void write(const Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<1>, int>& field, const std::string pathName);

  virtual void read(Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) const;
  virtual void read(Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) const;
  virtual void read(Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<1>, int>& field, const std::string pathName) const;
#endif

#ifdef SPHERAL2D
  virtual void write(const Field<Dim<2>, Dim<1>::Scalar>& field, const std::string pathName);
  virtual void write(const Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName);
  virtual void write(const Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<2>, int>& field, const std::string pathName);

  virtual void read(Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) const;
  virtual void read(Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) const;
  virtual void read(Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<2>, int>& field, const std::string pathName) const;
#endif

#ifdef SPHERAL3D
  virtual void write(const Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName);
  virtual void write(const Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName);
  virtual void write(const Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName);
  virtual void write(const Field<Dim<3>, int>& field, const std::string pathName);

  virtual void read(Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) const;
  virtual void read(Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) const;
  virtual void read(Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) const;
  virtual void read(Field<Dim<3>, int>& field, const std::string pathName) const;
#endif
  //******************************************************************************

  //------------------------------------------------------------------------------
  // We have to forward the templated write/read methods to the base class due to
  // function hiding.
  // Write/read FieldLists.
  template<typename Dimension, typename DataType>
  void write(const FieldList<Dimension, DataType>& fieldList, const std::string pathName) { FileIO::write(fieldList, pathName); }
  template<typename Dimension, typename DataType>
  void read(FieldList<Dimension, DataType>& fieldList, const std::string pathName) const { FileIO::read(fieldList, pathName); }

  // Write/read Fields of vectors.
  template<typename Dimension, typename DataType>
  void write(const Field<Dimension, std::vector<DataType> >& field, const std::string pathName) { FileIO::write(field, pathName); }
  template<typename Dimension, typename DataType>
  void read(Field<Dimension, std::vector<DataType> >& field, const std::string pathName) const {FileIO::read(field, pathName); }

  // Write/read a vector<DataType> if DataType is a primitive we already know about.
  template<typename DataType>
  void write(const std::vector<DataType>& x, const std::string pathName) { FileIO::write(x, pathName); }
  template<typename DataType>
  void read(std::vector<DataType>& x, const std::string pathName) const { FileIO::read(x, pathName); }
  //------------------------------------------------------------------------------

private:
  //--------------------------- Private Interface ---------------------------//
  // A pointer to the SiloFile file associated with this object.
  mutable DBfile* mFilePtr;

  // Don't allow assignment.
  SiloFileIO& operator=(const SiloFileIO& rhs);

  // Parse a path name and set the current working directory in a silo file.
  // The non-const version will create the directory in the silo file if need be,
  // while the const version will abort if the directory does not exist.
  std::string setDir(const std::string& pathName);
  std::string setDir(const std::string& pathName) const;

  // Helper methods to read/write sequential types.
  template<typename Container>
  void writeValueSequence(const Container& value, const std::string pathName, const int nvali);
  template<typename Container> 
  void readValueSequence(Container& value, const std::string pathName, const int nvali) const;
};

}

#else

// Forward declaration.
namespace Spheral {
  class SiloFileIO;
}

#endif
