//---------------------------------Spheral++----------------------------------//
// ArtificialConductionPolicy -- Override the default energy policy in the
// presence of artificial conduction.
//
// Created by CDR, 9/30/2014
//----------------------------------------------------------------------------//

#ifndef __ArtificialConductionPolicy_hh__
#define __ArtificialConductionPolicy_hh__

#include <string>
#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {
    
    // Forward declarations.
    template<typename Dimension> class State;
    template<typename Dimension> class StateDerivatives;
    template<typename Dimension> class FluidNodeList;
    template<typename Dimension, typename DataType> class Field;
    
    template<typename Dimension>
    class ArtificialConductionPolicy: public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
    public:
        //--------------------------- Public Interface ---------------------------//
        // Useful typedefs
        typedef typename Dimension::Scalar Scalar;
        typedef typename Dimension::SymTensor SymTensor;
        typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;
        typedef typename State<Dimension>::PolicyPointer PolicyPointer;
        
        // Constructors, destructor.
        ArtificialConductionPolicy(PolicyPointer& energyPolicy);
        virtual ~ArtificialConductionPolicy();
        
        // Overload the methods describing how to update Fields.
        virtual void update(const KeyType& key,
                            State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs,
                            const double multiplier,
                            const double t,
                            const double dt);
        
        virtual void updateAsIncrement(const KeyType& key,
                                       State<Dimension>& state,
                                       StateDerivatives<Dimension>& derivs,
                                       const double multiplier,
                                       const double t,
                                       const double dt);
        
        void conduct(const KeyType& key,
                     State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs,
                     const double multiplier,
                     const double t,
                     const double dt);
        
        // Equivalence.
        virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;
        
    private:
        //--------------------------- Private Interface ---------------------------//
        ArtificialConductionPolicy(const ArtificialConductionPolicy& rhs);
        ArtificialConductionPolicy& operator=(const ArtificialConductionPolicy& rhs);
        
        typename State<Dimension>::PolicyPointer mEnergyPolicy;
    };
    
}

#endif
