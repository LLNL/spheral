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
  virtual bool pathExists(const std::string path) const override;

  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const unsigned& value, const std::string path) override;
  virtual void write(const size_t& value, const std::string path) override;
  virtual void write(const int& value, const std::string path) override;
  virtual void write(const bool& value, const std::string path) override;
  virtual void write(const double& value, const std::string path) override;
  virtual void write(const std::string& value, const std::string path) override;
  virtual void write(const std::vector<int>& value, const std::string path) override;
  virtual void write(const std::vector<double>& value, const std::string path) override;
  virtual void write(const std::vector<std::string>& value, const std::string path) override;

  virtual void write(const Dim<1>::Vector& value, const std::string path) override;
  virtual void write(const Dim<1>::Tensor& value, const std::string path) override;
  virtual void write(const Dim<1>::SymTensor& value, const std::string path) override;
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string path) override;

  virtual void write(const Dim<2>::Vector& value, const std::string path) override;
  virtual void write(const Dim<2>::Tensor& value, const std::string path) override;
  virtual void write(const Dim<2>::SymTensor& value, const std::string path) override;
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string path) override;

  virtual void write(const Dim<3>::Vector& value, const std::string path) override;
  virtual void write(const Dim<3>::Tensor& value, const std::string path) override;
  virtual void write(const Dim<3>::SymTensor& value, const std::string path) override;
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string path) override;

  virtual void read(unsigned& value, const std::string path) const override;
  virtual void read(size_t& value, const std::string path) const override;
  virtual void read(int& value, const std::string path) const override;
  virtual void read(bool& value, const std::string path) const override;
  virtual void read(double& value, const std::string path) const override;
  virtual void read(std::string& value, const std::string path) const override;
  virtual void read(std::vector<int>& value, const std::string path) const override;
  virtual void read(std::vector<double>& value, const std::string path) const override;
  virtual void read(std::vector<std::string>& value, const std::string path) const override;

  virtual void read(Dim<1>::Vector& value, const std::string path) const override;
  virtual void read(Dim<1>::Tensor& value, const std::string path) const override;
  virtual void read(Dim<1>::SymTensor& value, const std::string path) const override;
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string path) const override;

  virtual void read(Dim<2>::Vector& value, const std::string path) const override;
  virtual void read(Dim<2>::Tensor& value, const std::string path) const override;
  virtual void read(Dim<2>::SymTensor& value, const std::string path) const override;
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string path) const override;

  virtual void read(Dim<3>::Vector& value, const std::string path) const override;
  virtual void read(Dim<3>::Tensor& value, const std::string path) const override;
  virtual void read(Dim<3>::SymTensor& value, const std::string path) const override;
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string path) const override;

  // //------------------------------------------------------------------------------
  // // We have to forward the templated write/read methods to the base class due to
  // // function hiding.
  // // Fields
  // template<typename Dimension, typename DataType> void write(const Field<Dimension, DataType>& value, const std::string path) const { FileIO::write(value, path); }
  // template<typename Dimension, typename DataType> void read(Field<Dimension, DataType>& value, const std::string path) { FileIO::read(value, path); }

  // // FieldLists
  // template<typename Dimension, typename DataType> void write(const FieldList<Dimension, DataType>& value, const std::string path) const { FileIO::write(value, path); }
  // template<typename Dimension, typename DataType> void read(FieldList<Dimension, DataType>& value, const std::string path) { FileIO::read(value, path); }

  // // Write/read a vector<Value> if Value is a primitive we already know about.
  // template<typename Value> void write(const std::vector<Value>& x, const std::string path) { FileIO::write(x, path); }
  // template<typename Value> void  read(std::vector<Value>& x, const std::string path) const { FileIO::read(x, path); }
  //------------------------------------------------------------------------------

  // Get and set the current precision for I/O.
  int precision() const;
  void setPrecision(int precision);

  // Test if we're currently prepared to write to the file.
  bool readyToWrite() const;

  // Test if we're currently prepared to read from the file.
  bool readyToRead() const;

  // Scan the file for the given path name.
  void findPathName(const std::string path) const;

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

  // Write generic DataTypes.
  template<typename DataType>
  void writeGenericType(const DataType& value,
                        const std::string path);

  // Read generic DataTypes.
  template<typename DataType>
  void readGenericType(DataType& value,
                       const std::string path) const;

  // Write a generic vector<T>.
  template<typename DataType>
  void writeGenericVector(const std::vector<DataType>& value,
                          const std::string path);

  // Read a generic vector<T>.
  template<typename DataType>
  void readGenericVector(std::vector<DataType>& value,
                         const std::string path) const;

};

}

#else

// Forward declaration.
namespace Spheral {
  class FlatFileIO;
}

#endif
