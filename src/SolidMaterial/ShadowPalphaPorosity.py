#-------------------------------------------------------------------------------
# ShadowPalphaPorosity
#
# Provides convenient constructors for the P-alpha porosity model
#-------------------------------------------------------------------------------
import types
from spheralDimensions import spheralDimensions
from CaptureStdout import helpString
from MaterialPropertiesLib import SpheralMaterialPropertiesLib
from SpheralCompiledPackages import *

dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic factory function
#-------------------------------------------------------------------------------
def _PalphaPorosityFactory(ndim):
    CXXPalphaPorosity = eval(f"PalphaPorosity{ndim}d")

    # The class we will return in place of the plain C++ one
    class PalphaPorosity(CXXPalphaPorosity):

        def __init__(self,
                     nodeList,
                     phi0,
                     Pe,
                     Pt,
                     Ps,
                     alphat,
                     n1,
                     n2,
                     cS0,
                     c0,
                     rho0 = None):
            """Construct a P-alpha porosity model.  Valid arguments are

  * nodeList: the (fluid or solid) NodeList whose nodes this porosity model will apply to
  *     phi0: (scalar or Field) initial porosity
  *       Pe: (scalar) Elastic pressure threshold
  *       Pt: (scalar) Transition pressure (Pe <= Pt)
  *       Ps: (scalar) Solid transition pressure (from fit, Pt <= Ps)
  *   alphat: (scalar) Transition distension
  *       n1: (scalar) Fitted exponent for plastic distention evolution
  *       n2: (scalar) Fitted exponent for plastic distention evolution
  *      cS0: (scalar) Reference sound speed at full density
  *       c0: (scalar or Field) Reference sound speed at initial porosity
  *     rho0: (scalar or None) Reference (solid) density.  If None, then reference value looked up from equation of state
"""

            if rho0 is None:
                try:
                    rho0 = nodeList.equationOfState.referenceDensity
                except:
                    raise RuntimeError("Unable to extract reference density for PalphaPorosity")
            CXXPalphaPorosity.__init__(self,
                                       nodeList = nodeList,
                                       phi0 = phi0,
                                       Pe = Pe,
                                       Pt = Pt,
                                       Ps = Ps,
                                       alphat = alphat,
                                       n1 = n1,
                                       n2 = n2,
                                       cS0 = cS0,
                                       c0 = c0,
                                       rho0 = rho0)
            return
    
    return PalphaPorosity

#-------------------------------------------------------------------------------
# Create the dimension specific Tillotson factories.  These are the ones
# you actually use.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec(f"PalphaPorosity{ndim}d = _PalphaPorosityFactor({ndim})")
