//
//  HelmholtzEquationOfState.cc
//  
//
//  Created by Raskin, Cody Dantes on 8/28/14.
//
//

#include "HelmholtzEquationOfState.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
    namespace Material {
        
        using FieldSpace::Field;
        

        template<typename Dimension>
        HelmholtzEquationOfState<Dimension>::
        HelmholtzEquationOfState(const PhysicalConstants& constants,
                                 const double minimumPressure,
                                 const double maximumPressure,
                                 const MaterialPressureMinType minPressureType,
                                 const Scalar abar0,
                                 const Scalar zbar0):
        EquationOfState<Dimension>(constants, minimumPressure, maximumPressure, minPressureType),
        mzbar(nullptr),
        mabar(nullptr),
        mabar0(abar0),
        mzbar0(zbar0)
        {
            
        }
        
        //------------------------------------------------------------------------------
        // Destructor.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        HelmholtzEquationOfState<Dimension>::
        ~HelmholtzEquationOfState() {
        }
        
        
        /* check mzbar != nullptr?? */
        
        
        //------------------------------------------------------------------------------
        // Set the pressure.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setPressure(Field<Dimension, Scalar>& Pressure,
                    const Field<Dimension, Scalar>& massDensity,
                    const Field<Dimension, Scalar>& specificThermalEnergy) const {
            CHECK(valid());
            for (size_t i = 0; i != massDensity.numElements(); ++i) {
                Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
            }
        }
        
        //------------------------------------------------------------------------------
        // Set the temperature.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setTemperature(Field<Dimension, Scalar>& temperature,
                       const Field<Dimension, Scalar>& massDensity,
                       const Field<Dimension, Scalar>& specificThermalEnergy) const {
            CHECK(valid());
            for (size_t i = 0; i != massDensity.numElements(); ++i) {
                temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
            }
        }
        
        //------------------------------------------------------------------------------
        // Set the specific thermal energy.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                 const Field<Dimension, Scalar>& massDensity,
                                 const Field<Dimension, Scalar>& temperature) const {
            CHECK(valid());
            for (size_t i = 0; i != massDensity.numElements(); ++i) {
                specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
            }
        }
        
        //------------------------------------------------------------------------------
        // Set the specific heat.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                        const Field<Dimension, Scalar>& massDensity,
                        const Field<Dimension, Scalar>& temperature) const {
            CHECK(valid());
            const double kB = mConstants.kB();
            const double mp = mConstants.protonMass();
            const double Cv = kB/(mGamma1*mMolecularWeight*mp);
            specificHeat = Cv;
        }
        
        //------------------------------------------------------------------------------
        // Set the sound speed.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                      const Field<Dimension, Scalar>& massDensity,
                      const Field<Dimension, Scalar>& specificThermalEnergy) const {
            CHECK(valid());
            for (size_t i = 0; i != massDensity.numElements(); ++i) {
                soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
            }
        }
        
        //------------------------------------------------------------------------------
        // Set gamma (ratio of specific heats).
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setGammaField(Field<Dimension, Scalar>& gamma,
                      const Field<Dimension, Scalar>& massDensity,
                      const Field<Dimension, Scalar>& specificThermalEnergy) const {
            CHECK(valid());
            gamma = mGamma;
        }
        
        //------------------------------------------------------------------------------
        // Set the bulk modulus (rho DP/Drho).  This is just the pressure for a
        // Helmholtz gas.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        void
        HelmholtzEquationOfState<Dimension>::
        setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                       const Field<Dimension, Scalar>& massDensity,
                       const Field<Dimension, Scalar>& specificThermalEnergy) const {
            CHECK(valid());
            setPressure(bulkModulus, massDensity, specificThermalEnergy);
            bulkModulus += mExternalPressure;
        }
        
        //------------------------------------------------------------------------------
        // Calculate an individual pressure.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        pressure(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const {
            CHECK(valid());
            return this->applyPressureLimits(mHelmholtzConstant*pow(massDensity, mGamma) - mExternalPressure);
        }
        
        //------------------------------------------------------------------------------
        // Calculate an individual temperature.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        temperature(const Scalar massDensity,
                    const Scalar specificThermalEnergy) const {
            CHECK(valid());
            const double kB = mConstants.kB();
            const double mp = mConstants.protonMass();
            return mGamma1*mMolecularWeight*mp/kB*specificThermalEnergy;
        }
        
        //------------------------------------------------------------------------------
        // Calculate an individual specific thermal energy.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        specificThermalEnergy(const Scalar massDensity,
                              const Scalar temperature) const {
            CHECK(valid());
            const double kB = mConstants.kB();
            return kB/(mGamma1*mMolecularWeight)*temperature;
        }
        
        //------------------------------------------------------------------------------
        // Calculate an individual specific heat.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        specificHeat(const Scalar massDensity,
                     const Scalar temperature) const {
            CHECK(valid());
            const double kB = mConstants.kB();
            const double mp = mConstants.protonMass();
            return kB/(mGamma1*mMolecularWeight*mp);
        }
        
        //------------------------------------------------------------------------------
        // Calculate an individual sound speed.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        soundSpeed(const Scalar massDensity,
                   const Scalar specificThermalEnergy) const {
            CHECK(valid());
            const double c2 = mHelmholtzConstant*pow(massDensity, mGamma1);
            CHECK(c2 >= 0.0);
            return sqrt(c2);
        }
        
        //------------------------------------------------------------------------------
        // Get gamma.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        gamma(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
            return mGamma;
        }
        
        //------------------------------------------------------------------------------
        // Calculate an individual bulk modulus.  
        // This is just the pressure for a Helmholtz gas.
        //------------------------------------------------------------------------------
        template<typename Dimension>
        typename Dimension::Scalar
        HelmholtzEquationOfState<Dimension>::
        bulkModulus(const Scalar massDensity,
                    const Scalar specificThermalEnergy) const {
            CHECK(valid());
            return pressure(massDensity, specificThermalEnergy) + mExternalPressure;
        }
        
        
        /* ACCESSORS */
        
        
        //------------------------------------------------------------------------------
        // Access abar
        //------------------------------------------------------------------------------
        template<typename Dimension>
        const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
        HelmholtzEquationOfState<Dimension>::
        abar() const {
            return mabar;
        }
        
        //------------------------------------------------------------------------------
        // Access zbar
        //------------------------------------------------------------------------------
        template<typename Dimension>
        const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
        HelmholtzEquationOfState<Dimension>::
        zbar() const {
            return mzbar;
        }
        
    }
}

