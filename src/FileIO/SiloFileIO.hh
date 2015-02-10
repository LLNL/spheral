//---------------------------------Spheral++----------------------------------//
// SiloFileIO -- Provide the interface to silo file objects.
//
// Created by JMO, Sat Feb  7 23:06:03 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_SiloFileIO__
#define __Spheral_SiloFileIO__

#include <vector>
#include <string>

#include "FileIO.hh"

extern "C" {
#include "silo.h"
}

namespace Spheral {
namespace FileIOSpace {

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
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<1>, int>& field, const std::string pathName);

  virtual void write(const FieldSpace::Field<Dim<2>, Dim<1>::Scalar>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<2>, int>& field, const std::string pathName);

  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName);
  virtual void write(const FieldSpace::Field<Dim<3>, int>& field, const std::string pathName);

  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<1>, int>& field, const std::string pathName) const;

  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<2>, int>& field, const std::string pathName) const;

  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) const;
  virtual void read(FieldSpace::Field<Dim<3>, int>& field, const std::string pathName) const;
  //******************************************************************************

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

};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace FileIOSpace {
    class SiloFileIO;
  }
}

#endif
