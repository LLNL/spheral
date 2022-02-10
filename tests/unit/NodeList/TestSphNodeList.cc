#include <iostream>
#include "NodeList/SphNodeList.hh"

void main () {

  typedef Dimension<3>::Scalar Scalar;
  typedef Dimension<3>::Vector Vector;
  typedef Dimension<3>::Tensor Tensor;
  typedef Dimension<3>::SymTensor SymTensor;

  // Create an SPH NodeList

  typedef Dimension<3>::Scalar Scalar;
  typedef Dimension<3>::Tensor Tensor;

  int nNodes = 10;
  cerr << "Creating an SPH node list of length " << nNodes << endl;

  SphNodeList< Dimension<3> > sphNodes(10);

  // Set the mass and velocity to known values.
  Scalar mass0 = 10.0;
  Vector velocity0(1.0, 2.0, 3.0);
  Field<Dimension<3>, Scalar> mass = sphNodes.mass();
  Field<Dimension<3>, Vector> velocity = sphNodes.velocity();
  for (int i = 0; i < mass.size(); ++i) {
    mass(i) = mass0;
    velocity(i) = velocity0;
  }

  cerr << "sphNodes.numNodes() = " << sphNodes.numNodes() << endl;
  cerr << "sphNodes.mass() = " << mass << endl;
  cerr << "sphNodes.velocity() = " << velocity << endl;

  Field<Dimension<3>, Vector> pmom = sphNodes.linearMomentum();
  cerr << "pmom = mass*velocity = " << pmom << endl;
  pmom = 10.0*pmom;
  cerr << "pmom = 10.0*pmom = " << pmom << endl;

  Field<Dimension<3>, Tensor> tensorField(sphNodes,
                                          Tensor(1, 2, 3,
                                                 4, 5, 6,
                                                 7, 8, 9));
  cerr << "tensorField = " << tensorField << endl;
  tensorField *= 10.0;
  cerr << "tensorField *= 10.0 => " << tensorField << endl;

  Field<Dimension<3>, Scalar> aScalarField(sphNodes);
  Field<Dimension<3>, Tensor> aTensorField(sphNodes);

  cerr << "aScalarField = " << aScalarField << endl;
  cerr << "aTensorField = " << aTensorField << endl;
  cerr << "aScalarField*aTensorField = " << aScalarField*aTensorField << endl;

}
