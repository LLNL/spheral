"""
Spheral DataBase module.

Provides higher level interface for State and DataBase objects.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# More preamble
#-------------------------------------------------------------------------------
preamble += """
PYBIND11_MAKE_OPAQUE(std::vector<FluidNodeList<Dim<1>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FluidNodeList<Dim<2>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FluidNodeList<Dim<3>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<NodeList<Dim<1>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<NodeList<Dim<2>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<NodeList<Dim<3>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<SolidNodeList<Dim<1>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<SolidNodeList<Dim<2>>*>)
PYBIND11_MAKE_OPAQUE(std::vector<SolidNodeList<Dim<3>>*>)

PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  int>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  int>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  int>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  int>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  int>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  int>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  int>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  int>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  int>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  int>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  int>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  int>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  double>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  double>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  double>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  double>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  double>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  double>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  double>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  double>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  double>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  double>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  double>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  double>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::SymTensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::SymTensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::SymTensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::SymTensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::SymTensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::SymTensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::Tensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::Tensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::Tensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::Tensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::Tensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::Tensor>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::Vector>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::Vector>*>)
PYBIND11_MAKE_OPAQUE(std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::Vector>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::Vector>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::Vector>*>)
PYBIND11_MAKE_OPAQUE(std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::Vector>*>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<1> PYB11COMMA  int>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<2> PYB11COMMA  int>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<3> PYB11COMMA  int>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<1> PYB11COMMA  int>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<2> PYB11COMMA  int>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<3> PYB11COMMA  int>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<1> PYB11COMMA  double>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<2> PYB11COMMA  double>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<3> PYB11COMMA  double>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<1> PYB11COMMA  double>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<2> PYB11COMMA  double>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<3> PYB11COMMA  double>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::SymTensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::SymTensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::SymTensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::SymTensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::SymTensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::SymTensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::Tensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::Tensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::Tensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::Tensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::Tensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::Tensor>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<1> PYB11COMMA  Dim<1>::Vector>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<2> PYB11COMMA  Dim<2>::Vector>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Field<Dim<3> PYB11COMMA  Dim<3>::Vector>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<1> PYB11COMMA  Dim<1>::Vector>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<2> PYB11COMMA  Dim<2>::Vector>>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<FieldList<Dim<3> PYB11COMMA  Dim<3>::Vector>>>)
"""

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes += ['"DataBase/DataBase.hh"',
             '"DataBase/StateBase.hh"',
             '"DataBase/State.hh"',
             '"DataBase/StateDerivatives.hh"',
             '"Field/Field.hh"',
             '"Neighbor/ConnectivityMap.hh"',
             '"Physics/Physics.hh"',
             '"Utilities/DataTypeTraits.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from StateBase import *
from State import *
from StateDerivatives import *
from DataBase import *

for ndim in dims:
    exec('''
StateBase%(ndim)id = PYB11TemplateClass(StateBase, template_parameters="Dim<%(ndim)i>")
State%(ndim)id = PYB11TemplateClass(State, template_parameters="Dim<%(ndim)i>")
StateDerivatives%(ndim)id = PYB11TemplateClass(StateDerivatives, template_parameters="Dim<%(ndim)i>")
DataBase%(ndim)id = PYB11TemplateClass(DataBase, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})
