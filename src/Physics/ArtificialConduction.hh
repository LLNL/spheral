//---------------------------------Spheral++----------------------------------//
// ArtificialConduction -- Artificial smoothing of energy discontinuities
//
//
// Created by CDR, 9/24/2014
//----------------------------------------------------------------------------//

#ifndef ArtificialConduction_HH
#define ArtificialConduction_HH

#include "Physics.hh"

namespace Spheral {
    template<typename Dimension> class State;
    template<typename Dimension> class StateDerivatives;
    
    namespace PhysicsSpace {
        template<typename Dimension>
        class ArtificialConduction: public Physics<Dimension> {
        public:
            //--------------------------- Public Interface ---------------------------//
            typedef typename Dimension::Scalar Scalar;
            typedef typename Dimension::Vector Vector;
            typedef typename Dimension::Tensor Tensor;
            typedef typename Dimension::SymTensor SymTensor;
            
            // Constructors
            ArtificialConduction(const bool gradPMode, const Scalar alphaArCond);
            
            // Destructor
            virtual ~ArtificialConduction();
            
            // Provide default methods for creating and registering an energy derivative.

            
            // Do any required one-time initializations on problem start up.
            virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);
            
            // Accessor Fns
            const FieldSpace::FieldList<Dimension, Scalar>& DepsDt() const;
            const FieldSpace::FieldList<Dimension, Scalar>& vsig() const;
            bool gradPMode() const;
            void gradPMode(bool val);
            
            // DepsDt Iterator Fn
            void computeConduction(const DataBaseSpace::DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state);
            
        private:
            //--------------------------- Private Interface ---------------------------//
            // Our derivative field(s).
            FieldSpace::FieldList<Dimension, Scalar> mDepsDt;
            FieldSpace::FieldList<Dimension, Scalar> mVsig;
            bool mGradPMode;
            Scalar mAlphaArCond;
            

        };
    }
}

#else

// Fwd Declaration
namespace Spheral{
    namespace PhysicsSpace {
        template<typename Dimension> class ArtificialConduction;
    }
}

#endif
