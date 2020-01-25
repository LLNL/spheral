//---------------------------------Spheral++----------------------------------//
// FVCRKHydroBase -- The FVCRKH/AFVCRKH hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Aug  9 14:30:26 PDT 2019
//----------------------------------------------------------------------------//
#ifndef __Spheral_FVCRKHydroBase_hh__
#define __Spheral_FVCRKHydroBase_hh__

#include "CRKSPH/CRKSPHHydroBase.hh"

namespace Spheral {

template<typename Dimension>
class FVCRKHydroBase: public CRKSPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  FVCRKHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  ArtificialViscosity<Dimension>& Q,
                  const TableKernel<Dimension>& W,
                  const TableKernel<Dimension>& WPi,
                  const double filter,
                  const double cfl,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool XSPH,
                  const MassDensityType densityUpdate,
                  const HEvolutionType HUpdate,
                  const RKOrder correctionOrder,
                  const double epsTensile,
                  const double nTensile,
                  const bool limitMultimaterialTopology);

  // Destructor.
  virtual ~FVCRKHydroBase();

  // // Tasks we do once on problem startup.
  // virtual
  // void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // // Register the state Hydro expects to use and evolve.
  // virtual 
  // void registerState(DataBase<Dimension>& dataBase,
  //                    State<Dimension>& state) override;

  // // Register the derivatives/change fields for updating state.
  // virtual
  // void registerDerivatives(DataBase<Dimension>& dataBase,
  //                          StateDerivatives<Dimension>& derivs) override;

  // // This method is called once at the beginning of a timestep, after all state registration.
  // virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
  //                                State<Dimension>& state,
  //                                StateDerivatives<Dimension>& derivs);

  // // Initialize the Hydro before we start a derivative evaluation.
  // virtual
  // void initialize(const Scalar time,
  //                 const Scalar dt,
  //                 const DataBase<Dimension>& dataBase,
  //                 State<Dimension>& state,
  //                 StateDerivatives<Dimension>& derivs) override;
                          
  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // // Finalize the derivatives.
  // virtual
  // void finalizeDerivatives(const Scalar time,
  //                          const Scalar dt,
  //                          const DataBase<Dimension>& dataBase,
  //                          const State<Dimension>& state,
  //                          StateDerivatives<Dimension>& derivs) const override;

  // // Finalize the hydro at the completion of an integration step.
  // virtual
  // void finalize(const Scalar time,
  //               const Scalar dt,
  //               DataBase<Dimension>& dataBase,
  //               State<Dimension>& state,
  //               StateDerivatives<Dimension>& derivs) override;
                  
  // // Apply boundary conditions to the physics specific fields.
  // virtual
  // void applyGhostBoundaries(State<Dimension>& state,
  //                           StateDerivatives<Dimension>& derivs) override;

  // // Enforce boundary conditions for the physics specific fields.
  // virtual
  // void enforceBoundaries(State<Dimension>& state,
  //                        StateDerivatives<Dimension>& derivs) override;

  // // // We need ghost connectivity to be computed.
  // // virtual bool requireGhostConnectivity() const override { return true; }

  // // Flag to choose whether we want to sum for density, or integrate
  // // the continuity equation.
  // MassDensityType densityUpdate() const;
  // void densityUpdate(MassDensityType type);

  // // Flag to select how we want to evolve the H tensor.
  // HEvolutionType HEvolution() const;
  // void HEvolution(HEvolutionType type);

  // // Flag to choose CRK Correction Order
  // RKOrder correctionOrder() const;
  // void correctionOrder(RKOrder order);

  // // Flag for the CRK volume weighting definition
  // CRKVolumeType volumeType() const;
  // void volumeType(CRKVolumeType x);

  // // Flag to determine if we're using the total energy conserving compatible energy
  // // evolution scheme.
  // bool compatibleEnergyEvolution() const;
  // void compatibleEnergyEvolution(bool val);

  // // Flag controlling if we evolve total or specific energy.
  // bool evolveTotalEnergy() const;
  // void evolveTotalEnergy(bool val);

  // // Flag to determine if we're using the XSPH algorithm.
  // bool XSPH() const;
  // void XSPH(bool val);

  // // Flag to determine if we cut multimaterial topology.
  // bool limitMultimaterialTopology() const;
  // void limitMultimaterialTopology(bool val);

  // // The object defining how we evolve smoothing scales.
  // const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // // Fraction of centroidal filtering to apply.
  // double filter() const;
  // void filter(double val);

  // // Parameters for the tensile correction force at small scales.
  // Scalar epsilonTensile() const;
  // void epsilonTensile(Scalar val);

  // Scalar nTensile() const;
  // void nTensile(Scalar val);
    
  // // We maintain a special boundary condition to handle void points.
  // const CRKSPHVoidBoundary<Dimension>& voidBoundary() const;

  // // The state field lists we're maintaining.
  // const FieldList<Dimension, int>&       timeStepMask() const;
  // const FieldList<Dimension, Scalar>&    pressure() const;
  // const FieldList<Dimension, Scalar>&    soundSpeed() const;
  // const FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  // const FieldList<Dimension, Scalar>&    entropy() const;
  // const FieldList<Dimension, SymTensor>& Hideal() const;
  // const FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  // const FieldList<Dimension, Scalar>&    effectiveViscousPressure() const;
  // const FieldList<Dimension, Scalar>&    viscousWork() const;
  // const FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  // const FieldList<Dimension, SymTensor>& massSecondMoment() const;
  // const FieldList<Dimension, Scalar>&    volume() const;
  // const FieldList<Dimension, Vector>&    massDensityGradient() const;
  // const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  // const FieldList<Dimension, Vector>&    DxDt() const;

  // const FieldList<Dimension, Vector>&    DvDt() const;
  // const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  // const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  // const FieldList<Dimension, SymTensor>& DHDt() const;
  // const FieldList<Dimension, Tensor>&    DvDx() const;
  // const FieldList<Dimension, Tensor>&    internalDvDx() const;
  // const FieldList<Dimension, std::vector<Vector> >& pairAccelerations() const;
  // const FieldList<Dimension, Vector>&    deltaCentroid() const;

  // const FieldList<Dimension, Scalar>&    A() const;
  // const FieldList<Dimension, Vector>&    B() const;
  // const FieldList<Dimension, Tensor>&    C() const;
  // const FieldList<Dimension, Vector>&    gradA() const;
  // const FieldList<Dimension, Tensor>&    gradB() const;
  // const FieldList<Dimension, ThirdRankTensor>&    gradC() const;
    
  // const FieldList<Dimension, Scalar>&                m0() const;
  // const FieldList<Dimension, Vector>&                m1() const;
  // const FieldList<Dimension, SymTensor>&             m2() const;
  // const FieldList<Dimension, ThirdRankTensor>&       m3() const;
  // const FieldList<Dimension, FourthRankTensor>&      m4() const;
  // const FieldList<Dimension, Vector>&                gradm0() const;
  // const FieldList<Dimension, Tensor>&                gradm1() const;
  // const FieldList<Dimension, ThirdRankTensor> &      gradm2() const;
  // const FieldList<Dimension, FourthRankTensor>&      gradm3() const;
  // const FieldList<Dimension, FifthRankTensor>&       gradm4() const;

  // const FieldList<Dimension, int>&       surfacePoint() const;
  // const FieldList<Dimension, std::vector<Vector>>& etaVoidPoints() const;
  // const FieldList<Dimension, FacetedVolume>& cells() const;
  // const FieldList<Dimension, std::vector<std::tuple<int, int, int>>>& cellFaceFlags() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "FVCRKHydroBase"; }
  // virtual void dumpState(FileIO& file, const std::string& pathName) const;
  // virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  FVCRKHydroBase();
  FVCRKHydroBase(const FVCRKHydroBase&);
  FVCRKHydroBase& operator=(const FVCRKHydroBase&);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class FVCRKHydroBase;
}

#endif
