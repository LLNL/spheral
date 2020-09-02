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
    
    //------------------------------------------------------------------------------
    // Constructor.
    //------------------------------------------------------------------------------
    template<typename Dimension>
    ArtificialConductionPolicy<Dimension>::
    ArtificialConductionPolicy(typename State<Dimension>::PolicyPointer& energyPolicy):
    FieldListUpdatePolicyBase<Dimension, Scalar>(),
    mEnergyPolicy(energyPolicy){
        const std::vector<std::string>& dependencies = energyPolicy->dependencies();
        for (unsigned i = 0; i<dependencies.size(); ++i)
            this->addDependency(dependencies[i]);
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
        mEnergyPolicy->update(key, state, derivs, multiplier, t, dt);
        conduct(key,state,derivs,multiplier,t,dt);
    }
    
    template<typename Dimension>
    void
    ArtificialConductionPolicy<Dimension>::
    updateAsIncrement(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) {
        mEnergyPolicy->updateAsIncrement(key,state,derivs,multiplier,t,dt);
        conduct(key,state,derivs,multiplier,t,dt);
    }
    
    template<typename Dimension>
    void
    ArtificialConductionPolicy<Dimension>::
    conduct(const KeyType& key,
            State<Dimension>& state,
            StateDerivatives<Dimension>& derivs,
            const double multiplier,
            const double /*t*/,
            const double /*dt*/) {
        KeyType fieldKey, nodeListKey;
        StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
        REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and
                nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
        
        FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
        const FieldList<Dimension, Scalar> DepsDt = derivs.fields("Artificial Cond DepsDt", 0.0);
        CHECK(eps.size() == DepsDt.size());

        // Loop over the internal values of the field.
        const unsigned numNodeLists = eps.size();
        for (unsigned k = 0; k != numNodeLists; ++k) {
            const unsigned n = eps[k]->numInternalElements();
            for (unsigned i = 0; i != n; ++i) {
                eps(k, i) += DepsDt(k,i) * multiplier;
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

