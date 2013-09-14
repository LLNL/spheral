//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//
// Created by JMO, Sun May 21 21:16:43 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosity__
#define __Spheral_ArtificialViscosity__

#include <utility>

#include "Geometry/Dimension.hh"

#ifndef __GCCXML__
#include <vector>
#include "Field/FieldList.hh"
#include "DataOutput/registerWithRestart.hh"
#else
#include "fakestl.hh"
#endif

namespace Spheral {

  template<typename Dimension> class AllNodeIterator;
  template<typename Dimension> class InternalNodeIterator;
  template<typename Dimension> class GhostNodeIterator;
  template<typename Dimension> class MasterNodeIterator;
  template<typename Dimension> class CoarseNodeIterator;
  template<typename Dimension> class RefineNodeIterator;

  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;

  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
    class FileIO;
  }
  namespace NeighborSpace {
    template<typename Dimension> class ConnectivityMap;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
}

namespace Spheral {
namespace ArtificialViscositySpace {

template<typename Dimension>
class ArtificialViscosity {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::pair<double, std::string> TimeStepType;

  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  // Constructors.
  ArtificialViscosity(const Scalar Clinear,
                      const Scalar Cquadratic);

  // Destructor.
  virtual ~ArtificialViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual void initialize(const DataBaseSpace::DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
			  const Scalar time, 
			  const Scalar dt,
                          const KernelSpace::TableKernel<Dimension>& W);

  // Require all descendents to return the artificial viscous Pi = P/rho^2 as a tensor.
  // Scalar viscosities should just return a diagonal tensor with their value along the diagonal.
  virtual std::pair<Tensor, Tensor> Piij(const unsigned nodeListi, const unsigned i, 
                                         const unsigned nodeListj, const unsigned j,
                                         const Vector& xi,
                                         const Vector& etai,
                                         const Vector& vi,
                                         const Scalar rhoi,
                                         const Scalar csi,
                                         const SymTensor& Hi,
                                         const Vector& xj,
                                         const Vector& etaj,
                                         const Vector& vj,
                                         const Scalar rhoj,
                                         const Scalar csj,
                                         const SymTensor& Hj) const = 0;

  // Allow access to the linear and quadratic constants.
  Scalar Cl() const;
  void Cl(Scalar Cl);

  Scalar Cq() const;
  void Cq(Scalar Cq);

  // Toggle for the Balsara shear correction.
  bool balsaraShearCorrection() const;
  void balsaraShearCorrection(bool value);

  // Calculate the curl of the velocity given the stress tensor.
  Scalar curlVelocityMagnitude(const Tensor& DvDx) const;

  // Access the FieldList of Balsara shear multiplicative corrections.
  const FieldSpace::FieldList<Dimension, Scalar>& shearMultiplier() const;

  // Access the internally computed estimate of sigma:
  // sig^ab = \partial v^a / \partial x^b.
  const FieldSpace::FieldList<Dimension, Tensor>& sigma() const;

  // Access the internally computed estimate of the velocity gradient and
  // grad div velocity.
  const FieldSpace::FieldList<Dimension, Vector>& gradDivVelocity() const;

  // Switch to turn the del^2 v limiter on/off.
  bool limiter() const;
  void limiter(bool value);

  // Access the epsilon safety factor.
  Scalar epsilon2() const;
  void epsilon2(Scalar epsilon2);

  // Method to return the limiter magnitude for the given node.
  Tensor calculateLimiter(const Vector& vi,
                          const Vector& vj,
                          const Scalar ci,
                          const Scalar cj,
                          const Scalar hi,
                          const Scalar hj,
                          const int nodeListID,
                          const int nodeID) const;

  // Helper for the limiter, calculate the unit grad div v term for the given 
  // node.
  Vector shockDirection(const Scalar ci,
                        const Scalar hi,
                        const int nodeListID,
                        const int nodeID) const;

  // The negligible sound speed parameter for use in the limiter.
  Scalar negligibleSoundSpeed() const;
  void negligibleSoundSpeed(const Scalar val);

  // The multiplier for sound speed in the limiter.
  Scalar csMultiplier() const;
  void csMultiplier(const Scalar val);

  // The multiplier for energy in the limiter.
  Scalar energyMultiplier() const;
  void energyMultiplier(const Scalar val);

  // Helper method to calculate Del cross V from the given sigma tensor.
  Scalar computeDelCrossVMagnitude(const Tensor& sigma) const;

  // Helper method to calculate the weighting based on the given position
  // for use in the sigma calculation.
  Vector sigmaWeighting(const Vector& r) const;

  // Figure out the total stress-strain tensor for a given node pair based on 
  // the stored average value and the given (position, velocity) pair.
  Tensor sigmaij(const Vector& rji, 
                 const Vector& rjiUnit, 
                 const Vector& vji, 
                 const Scalar& hi2,
                 const int nodeListID,
                 const int nodeID) const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "ArtificialViscosity"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Allow descendents to request that sigma be calculated.
  bool calculateSigma() const;
  void calculateSigma(bool value);

  Scalar mClinear;
  Scalar mCquadratic;

#ifndef __GCCXML__
  // Parameters for the Balsara shear correction.
  bool mBalsaraShearCorrection;
  FieldSpace::FieldList<Dimension, Scalar> mShearMultiplier;

  // Parameters for the Q limiter.
  bool mCalculateSigma;
  bool mLimiterSwitch;
  Scalar mEpsilon2;
  Scalar mNegligibleSoundSpeed;
  Scalar mCsMultiplier;
  Scalar mEnergyMultiplier;
  FieldSpace::FieldList<Dimension, Tensor> mSigma;
  FieldSpace::FieldList<Dimension, Vector> mGradDivVelocity;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif

  // Protected methods.
  virtual void calculateSigmaAndGradDivV(const DataBaseSpace::DataBase<Dimension>& dataBase,
                                         const State<Dimension>& state,
                                         const StateDerivatives<Dimension>& derivs,
                                         const KernelSpace::TableKernel<Dimension>& W,
                                         ConstBoundaryIterator boundaryBegin,
                                         ConstBoundaryIterator boundaryEnd);

private:
  //--------------------------- Private Interface ---------------------------//
  ArtificialViscosity();
  ArtificialViscosity(const ArtificialViscosity&);
  ArtificialViscosity& operator=(const ArtificialViscosity&) const;
};

}
}

#ifndef __GCCXML__
#include "ArtificialViscosityInline.hh"
#endif

#else

namespace Spheral {
  namespace ArtificialViscositySpace {
    // Forward declaration.
    template<typename Dimension> class ArtificialViscosity;
  }
}

#endif
