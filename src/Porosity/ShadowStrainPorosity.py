#-------------------------------------------------------------------------------
# ShadowStrainPorosity
#
# Provides convenient constructors for the P-alpha porosity model
#-------------------------------------------------------------------------------
import types
import scipy.integrate
from spheralDimensions import spheralDimensions
from CaptureStdout import helpString
from MaterialPropertiesLib import SpheralMaterialPropertiesLib
from PorousEquationOfState import *
from PorousStrengthModel import *
from SpheralCompiledPackages import *

dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Factory function to generate our Python enhanced StrainPorosity model.
#-------------------------------------------------------------------------------
def _StrainPorosityFactory(ndim):
    CXXStrainPorosity = eval(f"StrainPorosity{ndim}d")
    EquationOfState = eval(f"EquationOfState{ndim}d")
    StrengthModel = eval(f"StrengthModel{ndim}d")
    ScalarField = eval(f"ScalarField{ndim}d")

    # The class we will return in place of the plain C++ one
    class StrainPorosity(CXXStrainPorosity):

        #-----------------------------------------------------------------------
        # Constructor
        #-----------------------------------------------------------------------
        def __init__(self, *args, **kwargs):
            """Construct a strain-alpha porosity model.  Valid arguments are

  * nodeList: the solid NodeList whose nodes this porosity model will apply to
  *             phi0: (scalar or Field) initial porosity
  *             epsE: (scalar) Elastic-plastic transition strain
  *             epsX: (scalar) Threshold strain between compaction regimes
  *            kappa: (scalar) Compaction rate
  *          gammaS0: (scalar) Reference gamma at full density
  *              cS0: (scalar) Reference sound speed at full density
  *               c0: (scalar or Field) Reference sound speed at initial porosity
  *            rhoS0: (scalar or None) Reference (solid) density.  If None, then reference value looked up from equation of state
  * jutziStateUpdate: (bool, default True) Update state (deviatoric stress and damage) as described in Jutzi 2008
"""

            # Check for the deprecated porous EOS and strength models
            if len(args) > 2:
                if isinstance(args[0], EquationOfState):
                    print("DEPRECATION WARNING: StrainPorosity no longer requires a PorousEquationOfState -- ignoring")
                    args = args[1:]
                if isinstance(args[0], StrengthModel):
                    print("DEPRECATION WARNING: StrainPorosity no longer requires a PorousStrengthModel -- ignoring")
                    args = args[1:]

            # The valid arguments
            validArgs = ["nodeList",
                         "phi0",
                         "epsE",
                         "epsX",
                         "kappa",
                         "gammaS0",
                         "cS0",
                         "c0",
                         "rhoS0",
                         "jutziStateUpdate"]
            
            # Add any args to our kwargs dictionary
            assert len(args) <= len(validArgs)
            for i, arg in enumerate(args):
                kwargs[validArgs[i]] = arg

            # See if the deprecated interface was used in the keywords
            deprecatedArgs = ["porousEOS",
                              "porousStrength"]
            for darg in deprecatedArgs:
                if darg in kwargs:
                    print("DEPRECATION WARNING: StrainPorosity no longer accepts ", darg, " -- ignoring")
                    del kwargs[darg]

            assert "nodeList" in kwargs
            nodeList = kwargs["nodeList"]

            # Did we get a solid mass density?
            if "rhoS0" not in kwargs or kwargs["rhoS0"] is None:
                try:
                    kwargs["rhoS0"] = nodeList.equationOfState().referenceDensity
                except:
                    raise RuntimeError("Unable to extract reference density for StrainPorosity")

            if not "jutziStateUpdate" in kwargs:
                kwargs["jutziStateUpdate"] = True

            # Now build the actual C++ object
            CXXStrainPorosity.__init__(self, **kwargs)
            return

    return StrainPorosity

#-------------------------------------------------------------------------------
# Create the dimension specific factories.  These are the ones you actually use.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec(f"StrainPorosity{ndim}d = _StrainPorosityFactory({ndim})")
