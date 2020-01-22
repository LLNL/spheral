"""
Spheral DataBase module.

Provides higher level interface for State and DataBase objects.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Opaque types
#-------------------------------------------------------------------------------
PYB11opaque += ["std::vector<FluidNodeList<Dim<1>>*>",
                "std::vector<FluidNodeList<Dim<2>>*>",
                "std::vector<FluidNodeList<Dim<3>>*>",
                "std::vector<NodeList<Dim<1>>*>",
                "std::vector<NodeList<Dim<2>>*>",
                "std::vector<NodeList<Dim<3>>*>",
                "std::vector<SolidNodeList<Dim<1>>*>",
                "std::vector<SolidNodeList<Dim<2>>*>",
                "std::vector<SolidNodeList<Dim<3>>*>",

                "std::vector<Field<Dim<1> , int>>",
                "std::vector<Field<Dim<2> , int>>",
                "std::vector<Field<Dim<3> , int>>",
                "std::vector<FieldList<Dim<1> , int>>",
                "std::vector<FieldList<Dim<2> , int>>",
                "std::vector<FieldList<Dim<3> , int>>",
                "std::vector<FieldList<Dim<1> , int>*>",
                "std::vector<FieldList<Dim<2> , int>*>",
                "std::vector<FieldList<Dim<3> , int>*>",
                "std::vector<Field<Dim<1> , int>*>",
                "std::vector<Field<Dim<2> , int>*>",
                "std::vector<Field<Dim<3> , int>*>",
                "std::vector<Field<Dim<1> , double>>",
                "std::vector<Field<Dim<2> , double>>",
                "std::vector<Field<Dim<3> , double>>",
                "std::vector<FieldList<Dim<1> , double>>",
                "std::vector<FieldList<Dim<2> , double>>",
                "std::vector<FieldList<Dim<3> , double>>",
                "std::vector<FieldList<Dim<1> , double>*>",
                "std::vector<FieldList<Dim<2> , double>*>",
                "std::vector<FieldList<Dim<3> , double>*>",
                "std::vector<Field<Dim<1> , double>*>",
                "std::vector<Field<Dim<2> , double>*>",
                "std::vector<Field<Dim<3> , double>*>",
                "std::vector<Field<Dim<1> , Dim<1>::SymTensor>>",
                "std::vector<Field<Dim<2> , Dim<2>::SymTensor>>",
                "std::vector<Field<Dim<3> , Dim<3>::SymTensor>>",
                "std::vector<FieldList<Dim<1> , Dim<1>::SymTensor>>",
                "std::vector<FieldList<Dim<2> , Dim<2>::SymTensor>>",
                "std::vector<FieldList<Dim<3> , Dim<3>::SymTensor>>",
                "std::vector<FieldList<Dim<1> , Dim<1>::SymTensor>*>",
                "std::vector<FieldList<Dim<2> , Dim<2>::SymTensor>*>",
                "std::vector<FieldList<Dim<3> , Dim<3>::SymTensor>*>",
                "std::vector<Field<Dim<1> , Dim<1>::SymTensor>*>",
                "std::vector<Field<Dim<2> , Dim<2>::SymTensor>*>",
                "std::vector<Field<Dim<3> , Dim<3>::SymTensor>*>",
                "std::vector<Field<Dim<1> , Dim<1>::Tensor>>",
                "std::vector<Field<Dim<2> , Dim<2>::Tensor>>",
                "std::vector<Field<Dim<3> , Dim<3>::Tensor>>",
                "std::vector<FieldList<Dim<1> , Dim<1>::Tensor>>",
                "std::vector<FieldList<Dim<2> , Dim<2>::Tensor>>",
                "std::vector<FieldList<Dim<3> , Dim<3>::Tensor>>",
                "std::vector<FieldList<Dim<1> , Dim<1>::Tensor>*>",
                "std::vector<FieldList<Dim<2> , Dim<2>::Tensor>*>",
                "std::vector<FieldList<Dim<3> , Dim<3>::Tensor>*>",
                "std::vector<Field<Dim<1> , Dim<1>::Tensor>*>",
                "std::vector<Field<Dim<2> , Dim<2>::Tensor>*>",
                "std::vector<Field<Dim<3> , Dim<3>::Tensor>*>",
                "std::vector<Field<Dim<1> , Dim<1>::Vector>>",
                "std::vector<Field<Dim<2> , Dim<2>::Vector>>",
                "std::vector<Field<Dim<3> , Dim<3>::Vector>>",
                "std::vector<FieldList<Dim<1> , Dim<1>::Vector>>",
                "std::vector<FieldList<Dim<2> , Dim<2>::Vector>>",
                "std::vector<FieldList<Dim<3> , Dim<3>::Vector>>",
                "std::vector<FieldList<Dim<1> , Dim<1>::Vector>*>",
                "std::vector<FieldList<Dim<2> , Dim<2>::Vector>*>",
                "std::vector<FieldList<Dim<3> , Dim<3>::Vector>*>",
                "std::vector<Field<Dim<1> , Dim<1>::Vector>*>",
                "std::vector<Field<Dim<2> , Dim<2>::Vector>*>",
                "std::vector<Field<Dim<3> , Dim<3>::Vector>*>",
                "std::vector<std::vector<Field<Dim<1> , int>>>",
                "std::vector<std::vector<Field<Dim<2> , int>>>",
                "std::vector<std::vector<Field<Dim<3> , int>>>",
                "std::vector<std::vector<FieldList<Dim<1> , int>>>",
                "std::vector<std::vector<FieldList<Dim<2> , int>>>",
                "std::vector<std::vector<FieldList<Dim<3> , int>>>",
                "std::vector<std::vector<Field<Dim<1> , double>>>",
                "std::vector<std::vector<Field<Dim<2> , double>>>",
                "std::vector<std::vector<Field<Dim<3> , double>>>",
                "std::vector<std::vector<FieldList<Dim<1> , double>>>",
                "std::vector<std::vector<FieldList<Dim<2> , double>>>",
                "std::vector<std::vector<FieldList<Dim<3> , double>>>",
                "std::vector<std::vector<Field<Dim<1> , Dim<1>::SymTensor>>>",
                "std::vector<std::vector<Field<Dim<2> , Dim<2>::SymTensor>>>",
                "std::vector<std::vector<Field<Dim<3> , Dim<3>::SymTensor>>>",
                "std::vector<std::vector<FieldList<Dim<1> , Dim<1>::SymTensor>>>",
                "std::vector<std::vector<FieldList<Dim<2> , Dim<2>::SymTensor>>>",
                "std::vector<std::vector<FieldList<Dim<3> , Dim<3>::SymTensor>>>",
                "std::vector<std::vector<Field<Dim<1> , Dim<1>::Tensor>>>",
                "std::vector<std::vector<Field<Dim<2> , Dim<2>::Tensor>>>",
                "std::vector<std::vector<Field<Dim<3> , Dim<3>::Tensor>>>",
                "std::vector<std::vector<FieldList<Dim<1> , Dim<1>::Tensor>>>",
                "std::vector<std::vector<FieldList<Dim<2> , Dim<2>::Tensor>>>",
                "std::vector<std::vector<FieldList<Dim<3> , Dim<3>::Tensor>>>",
                "std::vector<std::vector<Field<Dim<1> , Dim<1>::Vector>>>",
                "std::vector<std::vector<Field<Dim<2> , Dim<2>::Vector>>>",
                "std::vector<std::vector<Field<Dim<3> , Dim<3>::Vector>>>",
                "std::vector<std::vector<FieldList<Dim<1> , Dim<1>::Vector>>>",
                "std::vector<std::vector<FieldList<Dim<2> , Dim<2>::Vector>>>",
                "std::vector<std::vector<FieldList<Dim<3> , Dim<3>::Vector>>>"]

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DataBase/DataBase.hh"',
                  '"DataBase/StateBase.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Field/Field.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"Physics/Physics.hh"',
                  '"Utilities/DataTypeTraits.hh"',
                  '"Geometry/CellFaceFlag.hh"',
                  '"RK/RKCoefficients.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

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
