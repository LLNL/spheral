//---------------------------------Spheral++----------------------------------//
// A simple form of the reducing artificial viscosity from Morris & Monaghan.
// Computes a correction scaling for the typical Monaghan & Gingold visc.
// References:
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//   Morris, J. P., & Monaghan, J. J. 1997, J. Comput. Phys., 136, 41
//
// Created by CDR, Aug 8th, 2014
//----------------------------------------------------------------------------//
#ifndef __MorrisMonaghanReducingViscosity__
#define __MorrisMonaghanReducingViscosity__

#include "ArtificialViscosity.hh"
#include "Physics/Physics.hh"

namespace Spheral {
namespace ArtificialViscositySpace {
    
template<typename Dimension>
    class MorrisMonaghanReducingViscosity: public ArtificialViscosity<Dimension> , public PhysicsSpace::Physics<Dimension>{
public:
    //--------------------------- Public Interface ---------------------------//
    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;
    typedef typename Dimension::SymTensor SymTensor;
    typedef typename PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;
    
    // Constructors & Destructors
    MorrisMonaghanReducingViscosity(const Scalar nh,
                                    const Scalar aStar,
                                    const Scalar Clinear,
                                    const Scalar Cquadratic,
                                    const bool linearInExpansion,
                                    const bool quadraticInExpansion);
    virtual ~MorrisMonaghanReducingViscosity();
    
    //............................................................................
    // Physics methods.
    // Increment the derivatives.
    virtual
    void evaluateDerivatives(const Scalar time,
                             const Scalar dt,
                             const DataBaseSpace::DataBase<Dimension>& dataBase,
                             const State<Dimension>& state,
                             StateDerivatives<Dimension>& derivatives) const;
    
    // Vote on a time step.
    virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase,
                            const State<Dimension>& state,
                            const StateDerivatives<Dimension>& derivs,
                            const Scalar currentTime) const;
    
    // Register our state.
    virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                               State<Dimension>& state);
    
    // Register the derivatives/change fields for updating state.
    virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                     StateDerivatives<Dimension>& derivs);
    
    // Do any required one-time initializations on problem start up.
    virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);
    //............................................................................
    
    // The required method to compute the artificial viscous P/rho^2.
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
                                           const SymTensor& Hj) const;
    
    // Access the switches for acting in expansion.
    bool linearInExpansion() const;
    void linearInExpansion(const bool x);
    
    bool quadraticInExpansion() const;
    void quadraticInExpansion(const bool x);
    
    // Restart methods.
    virtual std::string label() const { return "MorrisMonaghanReducingViscosity"; }
    
private:
    //--------------------------- Private Interface ---------------------------//
    bool mLinearInExpansion, mQuadraticInExpansion;
    
    MorrisMonaghanReducingViscosity();
    MorrisMonaghanReducingViscosity(const MorrisMonaghanReducingViscosity&);
    MorrisMonaghanReducingViscosity& operator=(const MorrisMonaghanReducingViscosity&) const;
    
    Scalar mnh,maStar;
    FieldSpace::FieldList<Dimension, Scalar> mrvAlpha, mDrvAlphaDt;
};
    
}
}

#else

namespace Spheral {
    namespace ArtificialViscositySpace {
        // Forward declaration.
        template<typename Dimension> class MorrisMonaghanReducingViscosity;
    }
}

#endif
