// BPL includes.
#include "boost/python.hpp"

namespace Spheral {

  // Declare the externally compiled extensions.
  void wrapInitializeTau();
  void wrapErff();
  void wrapNewtonRaphson();
  void wrapGlobalNodeIDs();
  void wrapRotationMatrix();
  void wrapIterateIdealH();
  void wrapCSPHUtilities();
  void wrapVisitWriter();
  void wrapMortonOrderIndicies();
  void wrapPeanoHilbertOrderIndicies();
  void wrapNodeOrdering();
  void wrapPackElement();

  BOOST_PYTHON_MODULE(Utilities) {
    wrapInitializeTau();
    wrapErff();
    wrapNewtonRaphson();
    wrapGlobalNodeIDs();
    wrapRotationMatrix();
    wrapIterateIdealH();
    wrapCSPHUtilities();
    wrapVisitWriter();
    wrapMortonOrderIndicies();
    wrapPeanoHilbertOrderIndicies();
    wrapNodeOrdering();
    wrapPackElement();
  }

}
