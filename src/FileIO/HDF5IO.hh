//---------------------------------Spheral++----------------------------------//
// HDF5IO -- Provide the interface to HDF5 file objects.
//
// Created by JMO, Tue Jul 11 22:19:37 PDT 2000
//----------------------------------------------------------------------------//
#ifndef HDF5IO_HH
#define HDF5IO_HH

#include "FileIO.hh"
#include "H5Cpp.h"

#include <string>
#include <map>

template<typename Dimension, typename DataType> class Field;

namespace Spheral {

template<typename Dimension>
class HDF5IO: public FileIO<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  HDF5IO();
  HDF5IO(const string& fileName, AccessType access);
  HDF5IO(const string& fileName, int access);

    // Destructor.
  virtual ~HDF5IO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const string& fileName, AccessType access);
  virtual void close();

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const int& value, const string& pathName);
  virtual void write(const bool& value, const string& pathName);
  virtual void write(const Scalar& value, const string& pathName);
  virtual void write(const Vector& value, const string& pathName);
  virtual void write(const Tensor& value, const string& pathName);
  virtual void write(const SymTensor& value, const string& pathName);
  virtual void write(const string& value, const string& pathName);

  virtual void read(int& value, const string& pathName) const;
  virtual void read(bool& value, const string& pathName) const;
  virtual void read(Scalar& value, const string& pathName) const;
  virtual void read(Vector& value, const string& pathName) const;
  virtual void read(Tensor& value, const string& pathName) const;
  virtual void read(SymTensor& value, const string& pathName) const;
  virtual void read(string& value, const string& pathName) const;

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
  virtual void write(const Field<Dimension, Scalar>& field, const string& pathName);
  virtual void write(const Field<Dimension, Vector>& field, const string& pathName);
  virtual void write(const Field<Dimension, Tensor>& field, const string& pathName);
  virtual void write(const Field<Dimension, SymTensor>& field, const string& pathName);

  virtual void read(Field<Dimension, Scalar>& field, const string& pathName) const;
  virtual void read(Field<Dimension, Vector>& field, const string& pathName) const;
  virtual void read(Field<Dimension, Tensor>& field, const string& pathName) const;
  virtual void read(Field<Dimension, SymTensor>& field, const string& pathName) const;
  //******************************************************************************

  // Write generic DataTypes.
  template<typename DataType, typename H5DataType>
  void writeGenericType(const DataType& value,
			const string& pathName,
			const H5DataType& h5Type);

  // Write generic containers.
  template<typename IteratorType, typename H5DataType>
  void writeGenericContainer(const IteratorType& begin,
                             const IteratorType& end,
			     const string& pathName,
                             const H5DataType& h5Type);

  // Write generic Fields.
  template<typename DataType, typename H5DataType>
  void writeGenericField(const Field<Dimension, DataType>& field,
			 const string& pathName,
                         const H5DataType& h5Type);

  // Read generic DataTypes.
  template<typename DataType, typename H5DataType>
  void readGenericType(DataType& value,
		       const string& pathName,
		       const H5DataType& h5Type) const;

  // Read generic containers.
  template<typename IteratorType, typename H5DataType>
  void readGenericContainer(const IteratorType& begin,
			    const IteratorType& end,
			    const string& pathName,
			    const H5DataType& h5Type) const;

  // Read generic Fields.
  template<typename DataType, typename H5DataType>
  void readGenericField(Field<Dimension, DataType>& field,
			const string& pathName,
			const H5DataType& h5Type) const;

  // Open the given group name, creating any groups necessary for the given path.
  Group* getH5Group(const string& pathName);

  // Test if we're currently prepared to write to the file.
  bool readyToWrite() const;

  // Test if we're currently prepared to read from the file.
  bool readyToRead() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // A static map relating FileIO AccessTypes with HDF5 types.
  map<AccessType, int> mHDF5AccessTypes;

  // A pointer to the HDF5 file associated with this object.
  H5File* mFilePtr;

  // Provide a method to set up the map of Spheral::FileIO access types to HDF5
  // file types.  This should be a static const so I don't have to do this!
  void initializeAccessMap();

  // Don't allow assignment.
  HDF5IO& operator=(const HDF5IO& rhs);

};

}

#endif
