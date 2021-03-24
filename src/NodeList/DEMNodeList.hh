//---------------------------------Spheral++----------------------------------//
// DEMNodeList -- 
// 3/10/2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMNodeList__
#define __Spheral_DEMNodeList__

#include "NodeList.hh"

#include <float.h>
#include <string>

namespace Spheral {

//template<typename Dimension> class NodeIteratorBase;
//template<typename Dimension> class AllNodeIterator;
//template<typename Dimension> class InternalNodeIterator;
//template<typename Dimension> class GhostNodeIterator;
//template<typename Dimension> class MasterNodeIterator;
//template<typename Dimension> class CoarseNodeIterator;
//template<typename Dimension> class RefineNodeIterator;
//template<typename Dimension> class State;
//template<typename Dimension> class StateDerivatives;
//template<typename Dimension> class Neighbor;
//template<typename Dimension> class ConnectivityMap;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
//template<typename Dimension> class TableKernel;
//template<typename Dimension> class DataBase;
class FileIO;

template<typename Dimension>
class DEMNodeList: public NodeList<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  //typedef typename Dimension::Tensor Tensor;
  //typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  DEMNodeList(std::string name,
                const int numInternal,
                const int numGhost,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const int maxNumNeighbors);

  // Destructor
  virtual ~DEMNodeList();

  // Access the fluid state variables.
  Field<Dimension, Vector>& angularVelocity();
  const Field<Dimension, Vector>& angularVelocity() const;
  void angularVelocity(const Field<Dimension, Vector>& angularVelocity);

  // These are quantities which are not stored, but can be computed.
//   virtual void volume(Field<Dimension, Scalar>& field) const;
//   virtual void linearMomentum(Field<Dimension, Vector>& field) const;
//   virtual void totalEnergy(Field<Dimension, Scalar>& field) const;


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
  // Fields that define the fluid's current state.
  Field<Dimension, Vector> mAngularVelocity;
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

