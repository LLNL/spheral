//---------------------------------Spheral++----------------------------------//
// DEMNodeList -- 
// 3/10/2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMNodeList__
#define __Spheral_DEMNodeList__

#include "NodeList.hh"
#include "Geometry/Dimension.hh"

#include <float.h>
#include <string>

namespace Spheral {

template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

class FileIO;

template<typename Dimension>
class DEMNodeList: public NodeList<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  // Constructors.
  DEMNodeList(std::string name,
                const int numInternal,
                const int numGhost,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const Scalar neighborSearchBuffer,
                const int maxNumNeighbors);

  // Destructor
  virtual ~DEMNodeList();

  Field<Dimension, Scalar>& particleRadius();
  const Field<Dimension, Scalar>& particleRadius() const;
  void particleRadius(const Field<Dimension, Scalar>& particleRadius);

  Field<Dimension, int>& compositeParticleIndex();
  const Field<Dimension, int>& compositeParticleIndex() const;
  void compositeParticleIndex(const Field<Dimension, int>& compositeParticleIndex);

  Field<Dimension, int>& uniqueIndex();
  const Field<Dimension, int>& uniqueIndex() const;
  void uniqueIndex(const Field<Dimension, int>& uniqueIndex);

  Scalar neighborSearchBuffer() const;
  void   neighborSearchBuffer(Scalar x);

  void setHfieldFromParticleRadius(const int startUniqueIndex);

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "DEMNodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  Scalar mNeighborSearchBuffer;
  Field<Dimension, Scalar> mParticleRadius;
  Field<Dimension, int> mCompositeParticleIndex;
  Field<Dimension, int> mUniqueIndex;
#endif

  // No default constructor or copying.
  DEMNodeList();
  DEMNodeList(const DEMNodeList& nodes);
  DEMNodeList& operator=(const DEMNodeList& rhs);
};

}

#include "DEMNodeListInline.hh"

#else

// Forward declaration
namespace Spheral {
  template<typename Dimension> class DEMNodeList;
}

#endif

