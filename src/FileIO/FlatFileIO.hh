//---------------------------------Spheral++----------------------------------//
// FlatFileIO -- Provide the interface to FlatFile file objects.
//
// Created by JMO, Fri Apr 13 01:19:02 PDT 2001
//----------------------------------------------------------------------------//
#ifndef FlatFileIO_HH
#define FlatFileIO_HH

#include "FileIO.hh"

#include <vector>
#include <string>
#include <fstream>

namespace Spheral {

enum class FlatFileFormat {
  ascii = 0,
  binary = 1
};

class FlatFileIO: public FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  FlatFileIO();
  FlatFileIO(const std::string fileName, AccessType access, FlatFileFormat format=FlatFileFormat::ascii);
  FlatFileIO(const std::string fileName, int access, int format=0);

  // Destructor.
  virtual ~FlatFileIO();

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

  //------------------------------------------------------------------------------
  // We have to forward the templated write/read methods to the base class due to
  // function hiding.
  // Write/read a vector<Value> if Value is a primitive we already know about.
  template<typename Value> void write(const std::vector<Value>& x, const std::string pathName) { FileIO::write(x, pathName); }
  template<typename Value> void  read(std::vector<Value>& x, const std::string pathName) const { FileIO::read(x, pathName); }
  //------------------------------------------------------------------------------

  // Get and set the current precision for I/O.
  int precision() const;
  void setPrecision(int precision);

  // Write generic DataTypes.
  template<typename DataType>
  void writeGenericType(const DataType& value,
                        const std::string pathName);

  // Read generic DataTypes.
  template<typename DataType>
  void readGenericType(DataType& value,
                       const std::string pathName) const;

  // Write a generic vector<T>.
  template<typename DataType>
  void writeGenericVector(const std::vector<DataType>& value,
                          const std::string pathName);

  // Read a generic vector<T>.
  template<typename DataType>
  void readGenericVector(std::vector<DataType>& value,
                         const std::string pathName) const;

  // Test if we're currently prepared to write to the file.
  bool readyToWrite() const;

  // Test if we're currently prepared to read from the file.
  bool readyToRead() const;

  // Scan the file for the given path name.
  void findPathName(const std::string pathName) const;

  // Move the current pointer to the beginning of the file.
  void beginningOfFile() const;

  // Forward hidden base class methods
  using FileIO::write;
  using FileIO::read;

private:
  //--------------------------- Private Interface ---------------------------//
  // The precision we will write to the output files with.
  int mPrecision;

  // A pointer to the FlatFile file associated with this object.
  mutable std::fstream* mFilePtr;

  FlatFileFormat mFileFormat;

  // Don't allow assignment.
  FlatFileIO& operator=(const FlatFileIO& rhs);

};

}

#else

// Forward declaration.
namespace Spheral {
  class FlatFileIO;
}

#endif
