//---------------------------------Spheral++----------------------------------//
// SidreFileIO -- Provide the interface to sidre file objects.
//
// Created by Mikhail Zakharchanka, 11/4/2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_SidreFileIO__
#define __Spheral_SidreFileIO__

#include "FileIO.hh"

#include <vector>

namespace Spheral
{

class SidreFileIO: public FileIO
{
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  SidreFileIO();
  SidreFileIO(const std::string fileName, AccessType access);

  // Destructor.
  virtual ~SidreFileIO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const std::string fileName, AccessType access) override;
  virtual void close() override;

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // Check if the specified path is in the file.
  virtual bool pathExists(const std::string pathName) const override;

  // All FileIO objects had better be able to read and write the primitive
  // DataTypes.
  virtual void write(const unsigned& value, const std::string pathName) override;
  virtual void write(const size_t& value, const std::string pathName) override;
  virtual void write(const int& value, const std::string pathName) override;
  virtual void write(const bool& value, const std::string pathName) override;
  virtual void write(const double& value, const std::string pathName) override;
  virtual void write(const std::string& value, const std::string pathName) override;
  virtual void write(const std::vector<int>& value, const std::string pathName) override;
  virtual void write(const std::vector<double>& value, const std::string pathName) override;
  virtual void write(const std::vector<std::string>& value, const std::string pathName) override;

  virtual void write(const Dim<1>::Vector& value, const std::string pathName) override;
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName) override;
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName) override;
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName) override;

  virtual void write(const Dim<2>::Vector& value, const std::string pathName) override;
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName) override;
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName) override;
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName) override;

  virtual void write(const Dim<3>::Vector& value, const std::string pathName) override;
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName) override;
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName) override;
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName) override;

  virtual void read(unsigned& value, const std::string pathName) const override;
  virtual void read(size_t& value, const std::string pathName) const override;
  virtual void read(int& value, const std::string pathName) const override;
  virtual void read(bool& value, const std::string pathName) const override;
  virtual void read(double& value, const std::string pathName) const override;
  virtual void read(std::string& value, const std::string pathName) const override;
  virtual void read(std::vector<int>& value, const std::string pathName) const override;
  virtual void read(std::vector<double>& value, const std::string pathName) const override;
  virtual void read(std::vector<std::string>& value, const std::string pathName) const override;

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const override;
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const override;
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const override;
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const override;

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const override;
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const override;
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const override;
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const override;

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const override;
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const override;
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const override;
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const override;

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
#ifdef SPHERAL1D
  virtual void write(const Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<1>, int>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<1>, unsigned>& field, const std::string pathName) override;

  virtual void read(Field<Dim<1>, Dim<1>::Scalar>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<1>, Dim<1>::Vector>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<1>, Dim<1>::Tensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<1>, Dim<1>::SymTensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<1>, int>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<1>, unsigned>& field, const std::string pathName) const override;
#endif

#ifdef SPHERAL2D
  virtual void write(const Field<Dim<2>, Dim<1>::Scalar>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<2>, int>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<2>, unsigned>& field, const std::string pathName) override;

  virtual void read(Field<Dim<2>, Dim<2>::Scalar>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<2>, Dim<2>::Vector>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<2>, Dim<2>::Tensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<2>, Dim<2>::SymTensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<2>, int>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<2>, unsigned>& field, const std::string pathName) const override;
#endif

#ifdef SPHERAL3D
  virtual void write(const Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<3>, int>& field, const std::string pathName) override;
  virtual void write(const Field<Dim<3>, unsigned>& field, const std::string pathName) override;

  virtual void read(Field<Dim<3>, Dim<3>::Scalar>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<3>, Dim<3>::Vector>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<3>, Dim<3>::Tensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<3>, Dim<3>::SymTensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<3>, int>& field, const std::string pathName) const override;
  virtual void read(Field<Dim<3>, unsigned>& field, const std::string pathName) const override;
#endif
  //******************************************************************************

  // Forward hidden base class methods
  using FileIO::write;
  using FileIO::read;

private:
  //--------------------------- Private Interface ---------------------------//
  // A pointer to the root of the sidre datastore associated with this object.
  std::shared_ptr<axom::sidre::DataStore> mDataStorePtr;

  //save() function in sidre needs to have access to file name, also used for loadExternalData()
  std::string mFileName;

  // Don't allow assignment.
  SidreFileIO& operator=(const SidreFileIO& rhs);
};

}

#else

// Forward declaration.
namespace Spheral
{
  class SidreFileIO;
}

#endif
