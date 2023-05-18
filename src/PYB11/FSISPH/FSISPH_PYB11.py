"""
Spheral FSISPH module.

"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()


from SolidFSISPHHydroBase import *
from SlideSurface import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"FSISPH/SolidFSISPHHydroBase.hh"',
                  '"FSISPH/FSIFieldNames.hh"',
                  '"FSISPH/SlideSurface.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"']

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
InterfaceMethod = PYB11enum(("HLLCInterface", 
                             "ModulusInterface",
                             "NoInterface"), export_values = True)

KernelAveragingMethod = PYB11enum(("NeverAverageKernels", 
                                   "AlwaysAverageKernels",
                                   "AverageInterfaceKernels"), export_values = True)

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
SlideSurface%(ndim)id = PYB11TemplateClass(SlideSurface, template_parameters="%(Dimension)s")
SolidFSISPHHydroBase%(ndim)id = PYB11TemplateClass(SolidFSISPHHydroBase, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})


#-------------------------------------------------------------------------------
# expose our field names
#-------------------------------------------------------------------------------
class FSIFieldNames:
    rawPressure = PYB11readonly(static=True, returnpolicy="copy")
    pressureGradient = PYB11readonly(static=True, returnpolicy="copy")
    specificThermalEnergyGradient = PYB11readonly(static=True, returnpolicy="copy")
    interfaceNormals = PYB11readonly(static=True, returnpolicy="copy")
    interfaceFraction = PYB11readonly(static=True, returnpolicy="copy")
    interfaceSmoothness = PYB11readonly(static=True, returnpolicy="copy")
    smoothedInterfaceNormals = PYB11readonly(static=True, returnpolicy="copy")
    smoothnessNormalization = PYB11readonly(static=True, returnpolicy="copy")
