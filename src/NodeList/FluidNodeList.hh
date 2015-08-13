//---------------------------------Spheral++----------------------------------//
// FluidNodeList -- An abstract base class for the NodeLists to represent
//                  fluids.
//
// Created by JMO, Sun Sep 12 17:12:24 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_FluidNodeList__
#define __Spheral_FluidNodeList__

#include <float.h>
#include <string>

#include "NodeList.hh"

namespace Spheral {
  template<typename Dimension> class NodeIteratorBase;
  template<typename Dimension> class AllNodeIterator;
  template<typename Dimension> class InternalNodeIterator;
  template<typename Dimension> class GhostNodeIterator;
  template<typename Dimension> class MasterNodeIterator;
  template<typename Dimension> class CoarseNodeIterator;
  template<typename Dimension> class RefineNodeIterator;
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NeighborSpace {
    template<typename Dimension> class Neighbor;
    template<typename Dimension> class ConnectivityMap;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace Material {
    template<typename Dimension> class EquationOfState;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FileIOSpace {
    class FileIO;
  }
  namespace ArtificialViscositySpace {
    template<typename Dimension> class ArtificialViscosity;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {
namespace NodeSpace {

template<typename Dimension>
class FluidNodeList: public NodeList<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  FluidNodeList(std::string name,
                Material::EquationOfState<Dimension>& eos,
                const int numInternal,
                const int numGhost,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const int maxNumNeighbors,
                const Scalar rhoMin,
                const Scalar rhoMax);

  // Destructor
  virtual ~FluidNodeList();

  // Access the fluid state variables.
  FieldSpace::Field<Dimension, Scalar>& massDensity();
  FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy();

  const FieldSpace::Field<Dimension, Scalar>& massDensity() const;
  const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy() const;

  void massDensity(const FieldSpace::Field<Dimension, Scalar>& rho);
  void specificThermalEnergy(const FieldSpace::Field<Dimension, Scalar>& eps);

  // These are quantities which are not stored, but can be computed.
  virtual void pressure(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void temperature(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void soundSpeed(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void volume(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void linearMomentum(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void totalEnergy(FieldSpace::Field<Dimension, Scalar>& field) const;

  // Access the equation of state.
  const Material::EquationOfState<Dimension>& equationOfState() const;
  void equationOfState(const Material::EquationOfState<Dimension>& eos);

  // Optional bounding mass densities for use when time integrating
  // the density.
  Scalar rhoMin() const;
  Scalar rhoMax() const;

  void rhoMin(const Scalar x);
  void rhoMax(const Scalar x);

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "FluidNodeList"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Min/max mass densities.
  Scalar mRhoMin, mRhoMax;

#ifndef __GCCXML__
  // Fields that define the fluid's current state.
  FieldSpace::Field<Dimension, Scalar> mMassDensity;
  FieldSpace::Field<Dimension, Scalar> mSpecificThermalEnergy;

  // Equation of state.
  const Material::EquationOfState<Dimension>* mEosPtr;
#endif

  // No default constructor or copying.
  FluidNodeList();
  FluidNodeList(const FluidNodeList& nodes);
  FluidNodeList& operator=(const FluidNodeList& rhs);
};

}
}

#ifndef __GCCXML__
#include "FluidNodeListInline.hh"
#endif

#else

// Forward declaration
namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class FluidNodeList;
  }
}

#endif

