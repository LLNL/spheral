//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//

#ifndef ArtificialConduction_HH
#define ArtificialConduction_HH

#include "Physics/Physics.hh"
#include "RK/RKCorrectionParams.hh"

namespace Spheral {
    template<typename Dimension> class State;
    template<typename Dimension> class StateDerivatives;
    
    template<typename Dimension>
    class ArtificialConduction: public Physics<Dimension> {
    public:
        //--------------------------- Public Interface ---------------------------//
        typedef typename Dimension::Scalar Scalar;
        typedef typename Dimension::Vector Vector;
        typedef typename Dimension::Tensor Tensor;
        typedef typename Dimension::SymTensor SymTensor;
        typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
            
        typedef typename Physics<Dimension>::TimeStepType TimeStepType;
            
        // Constructors
        ArtificialConduction(const TableKernel<Dimension>& W,
                             const Scalar alphaArCond, const RKOrder ACcorrectionOrder = RKOrder::LinearOrder);
            
        // Destructor
        virtual ~ArtificialConduction();
            
        // Do any required one-time initializations on problem start up.
        virtual void initializeProblemStartup(DataBase<Dimension>& dataBase);
            
        // Register our state.
        virtual void registerState(DataBase<Dimension>& dataBase,
                                   State<Dimension>& state);
        //Allow access to the AC correction order.
        RKOrder ACcorrectionOrder() const;
        void ACcorrectionOrder(RKOrder order);
            
        // Provide default methods for registering and iterating derivatives.
        virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                         StateDerivatives<Dimension>& derivs);
        virtual
        void evaluateDerivatives(const Scalar time,
                                 const Scalar dt,
                                 const DataBase<Dimension>& dataBase,
                                 const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivatives) const;
            
            
        // Vote on a time step.
        virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                                const State<Dimension>& state,
                                const StateDerivatives<Dimension>& derivs,
                                const Scalar currentTime) const;
            
        virtual std::string label() const { return "Artificial Conduction"; }
            
        // Accessor Fns
            
    private:
        //--------------------------- Private Interface ---------------------------//
        const TableKernel<Dimension>& mKernel;
            
        // Our derivative field(s).
        FieldList<Dimension, Vector> mGradP;
        FieldList<Dimension, Scalar> mDepsDtArty;
        FieldList<Dimension, Scalar> mVsigMax;
        Scalar mAlphaArCond;
        RKOrder mACcorrectionOrder;

    };
}

#else

// Fwd Declaration
namespace Spheral{
    template<typename Dimension> class ArtificialConduction;
}

#endif
