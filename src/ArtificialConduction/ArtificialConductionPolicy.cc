//---------------------------------Spheral++----------------------------------//
// ArtificialConductionPolicy -- Override the default energy policy in the
// presence of artificial conduction.
//----------------------------------------------------------------------------//
#include "ArtificialConductionPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
    
    using NodeSpace::NodeList;
    using NodeSpace::FluidNodeList;
    using FieldSpace::FieldList;
    
    //------------------------------------------------------------------------------
    // Constructor.
    //------------------------------------------------------------------------------
    template<typename Dimension>
    ArtificialConductionPolicy<Dimension>::
    ArtificialConductionPolicy(State<Dimension>::PolicyPointer& energyPolicy):
    FieldListUpdatePolicyBase<Dimension>(),
    mEnergyPolicy(energyPolicy){
        
    }
    
    //------------------------------------------------------------------------------
    // Destructor.
    //------------------------------------------------------------------------------
    template<typename Dimension>
    ArtificialConductionPolicy<Dimension>::
    ~ArtificialConductionPolicy() {
    }
    
    //------------------------------------------------------------------------------
    // Update the field.
    //------------------------------------------------------------------------------
    template<typename Dimension>
    void
    ArtificialConductionPolicy<Dimension>::
    update(const KeyType& key,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs,
           const double multiplier,
           const double t,
           const double dt) {
        KeyType fieldKey, nodeListKey;
        StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
        REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and
                nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
        
        
        
        FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
        const FieldSpace::FieldList<Dimension, Value> DepsDt = derivs.fields("Artificial Cond. DepsDt", Scalar);
        CHECK(eps.size() == DepsDt.size());
        
        // Have the base class update the energy.
        mEnergyPolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);
        
        // Get the artycond depsdt from derivs

        
        // Loop over the internal values of the field.
        const unsigned numNodeLists = eps.size();
        for (unsigned k = 0; k != numNodeLists; ++k) {
            const unsigned n = eps[k]->numInternalElements();
            for (unsigned i = 0; i != n; ++i) {
                eps(k, i) += DepsDt(k,i) * dt;
            }
        }

    }
    
    //------------------------------------------------------------------------------
    // Equivalence operator.
    //------------------------------------------------------------------------------
    template<typename Dimension>
    bool
    ArtificialConductionPolicy<Dimension>::
    operator==(const UpdatePolicyBase<Dimension>& rhs) const {
        
        // We're only equal if the other guy is also an increment operator.
        const ArtificialConductionPolicy<Dimension>* rhsPtr = dynamic_cast<const ArtificialConductionPolicy<Dimension>*>(&rhs);
        if (rhsPtr == 0) {
            return false;
        } else {
            return true;
        }
    }
    
}

