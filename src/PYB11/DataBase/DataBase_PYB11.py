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
                  '"DataBase/UpdatePolicyBase.hh"',
                  '"DataBase/FieldUpdatePolicyBase.hh"',
                  '"DataBase/CopyState.hh"',
                  '"DataBase/IncrementState.hh"',
                  '"DataBase/IncrementBoundedState.hh"',
                  '"DataBase/ReplaceState.hh"',
                  '"DataBase/ReplaceBoundedState.hh"',
                  '"DataBase/FieldListUpdatePolicyBase.hh"',
                  '"DataBase/IncrementFieldList.hh"',
                  '"DataBase/IncrementBoundedFieldList.hh"',
                  '"DataBase/ReplaceFieldList.hh"',
                  '"DataBase/ReplaceBoundedFieldList.hh"',
                  '"DataBase/CompositeFieldListPolicy.hh"',
                  '"DataBase/CopyFieldList.hh"',
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
from UpdatePolicyBase import *
from FieldUpdatePolicyBase import *
from CopyState import *
from IncrementState import *
from IncrementBoundedState import *
from ReplaceState import *
from ReplaceBoundedState import *
from FieldListUpdatePolicyBase import *
from IncrementFieldList import *
from IncrementBoundedFieldList import *
from ReplaceFieldList import *
from ReplaceBoundedFieldList import *
from CompositeFieldListPolicy import *
from CopyFieldList import *

for ndim in dims:
    exec(f'''
StateBase{ndim}d = PYB11TemplateClass(StateBase, template_parameters="Dim<{ndim}>")
State{ndim}d = PYB11TemplateClass(State, template_parameters="Dim<{ndim}>")
StateDerivatives{ndim}d = PYB11TemplateClass(StateDerivatives, template_parameters="Dim<{ndim}>")
DataBase{ndim}d = PYB11TemplateClass(DataBase, template_parameters="Dim<{ndim}>")
UpdatePolicyBase{ndim}d = PYB11TemplateClass(UpdatePolicyBase, template_parameters="Dim<{ndim}>")
''')

    for (value, label) in (("Dim<%i>::Scalar" % ndim,        "Scalar"),
                          ):
        Dimension = f"Dim<{ndim}>"
        suffix = f"{ndim}d"
        exec(f'''
{label}FieldUpdatePolicyBase{suffix} = PYB11TemplateClass(FieldUpdatePolicyBase, template_parameters = ("{Dimension}", "{value}"))
{label}CopyState{suffix} = PYB11TemplateClass(CopyState, template_parameters = ("{Dimension}", "{value}"))
{label}IncrementState{suffix} = PYB11TemplateClass(IncrementState, template_parameters = ("{Dimension}", "{value}"))
{label}IncrementBoundedState{suffix} = PYB11TemplateClass(IncrementBoundedState, template_parameters = ("{Dimension}", "{value}"))
{label}ReplaceState{suffix} = PYB11TemplateClass(ReplaceState, template_parameters = ("{Dimension}", "{value}"))
{label}ReplaceBoundedState{suffix} = PYB11TemplateClass(ReplaceBoundedState, template_parameters = ("{Dimension}", "{value}"))
{label}FieldListUpdatePolicyBase{suffix} = PYB11TemplateClass(FieldListUpdatePolicyBase, template_parameters = ("{Dimension}", "{value}"))
{label}IncrementFieldList{suffix} = PYB11TemplateClass(IncrementFieldList, template_parameters = ("{Dimension}", "{value}"))
{label}IncrementBoundedFieldList{suffix} = PYB11TemplateClass(IncrementBoundedFieldList, template_parameters = ("{Dimension}", "{value}"))
{label}ReplaceFieldList{suffix} = PYB11TemplateClass(ReplaceFieldList, template_parameters = ("{Dimension}", "{value}"))
{label}ReplaceBoundedFieldList{suffix} = PYB11TemplateClass(ReplaceBoundedFieldList, template_parameters = ("{Dimension}", "{value}"))
{label}CompositeFieldListPolicy{suffix} = PYB11TemplateClass(CompositeFieldListPolicy, template_parameters = ("{Dimension}", "{value}"))
{label}CopyFieldList{suffix} = PYB11TemplateClass(CopyFieldList, template_parameters = ("{Dimension}", "{value}"))
''')
