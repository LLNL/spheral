//---------------------------------Spheral++----------------------------------//
// SolidNodeList -- A form of the FluidNodeList appropriate for use with 
// solid materials.
//
// Created by JMO, Tue Sep 7 22:44:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidNodeList_hh__
#define __Spheral_SolidNodeList_hh__

#include "NodeList/FluidNodeList.hh"

namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class StrengthModel;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class SolidNodeList: public NodeSpace::FluidNodeList<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  SolidNodeList(std::string name,
                Material::EquationOfState<Dimension>& eos,
                StrengthModel<Dimension>& strength,
                const int numInternal,
                const int numGhost,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const int maxNumNeighbors,
                const Scalar rhoMin,
                const Scalar rhoMax);

  // Destructor.
  virtual ~SolidNodeList();

  // Override the base method for calculating the sound speed.
  virtual void soundSpeed(FieldSpace::Field<Dimension, Scalar>& field) const;

  // Compute the bulk modulus, shear modulus, and yield strength.
  virtual void bulkModulus(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void shearModulus(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void yieldStrength(FieldSpace::Field<Dimension, Scalar>& field) const;

  // Access the stored deviatoric stress and plastic strain.
  FieldSpace::Field<Dimension, SymTensor>& deviatoricStress();
  FieldSpace::Field<Dimension, Scalar>& plasticStrain();
  FieldSpace::Field<Dimension, Scalar>& plasticStrainRate();

  const FieldSpace::Field<Dimension, SymTensor>& deviatoricStress() const;
  const FieldSpace::Field<Dimension, Scalar>& plasticStrain() const;
  const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate() const;

  // The tensor damage and it's gradient.
  FieldSpace::Field<Dimension, SymTensor>& damage();
  FieldSpace::Field<Dimension, SymTensor>& effectiveDamage();
  FieldSpace::Field<Dimension, Vector>& damageGradient();

  const FieldSpace::Field<Dimension, SymTensor>& damage() const;
  const FieldSpace::Field<Dimension, SymTensor>& effectiveDamage() const;
  const FieldSpace::Field<Dimension, Vector>& damageGradient() const;

  // The fragment ID field.
  FieldSpace::Field<Dimension, int>& fragmentIDs();
  const FieldSpace::Field<Dimension, int>& fragmentIDs() const;

  // The strength model this solid is using.
  const SolidMaterial::StrengthModel<Dimension>& strengthModel() const;

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "SolidNodeList"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  FieldSpace::Field<Dimension, SymTensor> mDeviatoricStress;
  FieldSpace::Field<Dimension, Scalar> mPlasticStrain;
  FieldSpace::Field<Dimension, Scalar> mPlasticStrainRate;
  FieldSpace::Field<Dimension, SymTensor> mDamage;
  FieldSpace::Field<Dimension, SymTensor> mEffectiveDamage;
  FieldSpace::Field<Dimension, Vector> mDamageGradient;
  FieldSpace::Field<Dimension, int> mFragmentIDs;
#endif

  // Pointer to the associated strength object.
  StrengthModel<Dimension>& mStrength;

  // No default constructor or copying.
  SolidNodeList();
  SolidNodeList(const SolidNodeList& solidNodeList);
};

}
}

#ifndef __GCCXML__
#include "SolidNodeListInline.hh"
#endif

#else

namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class SolidNodeList;
  }
}

#endif
