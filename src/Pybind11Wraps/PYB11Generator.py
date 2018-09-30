#-------------------------------------------------------------------------------
# PYB11Generator
#-------------------------------------------------------------------------------
import inspect
import sys
from PYB11utils import *
from PYB11Decorators import *
from PYB11STLmethods import *
from PYB11function import *
from PYB11class import *
from PYB11enum import *
from PYB11attr import *

#-------------------------------------------------------------------------------
# PYB11generateModule
#-------------------------------------------------------------------------------
def PYB11generateModule(modobj):
    name = modobj.__name__
    with open(name + ".cc", "w") as f:
        ss = f.write
        PYB11generateModuleStart(modobj, ss, name)

        # enums
        PYB11generateModuleEnums(modobj, ss)

        # classes
        PYB11generateModuleClasses(modobj, ss)

        # STL types
        PYB11generateModuleSTL(modobj, ss)

        # methods
        PYB11generateModuleFunctions(modobj, ss)

        # Attributes
        PYB11generateModuleAttrs(modobj, ss)

        # Closing
        ss("}\n")
        f.close()

    return

#-------------------------------------------------------------------------------
# PYB11generateModuleStart
#
# All the stuff up to the methods.
#-------------------------------------------------------------------------------
def PYB11generateModuleStart(modobj, ss, name):

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Module %(name)s
//------------------------------------------------------------------------------
// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"
#include "pybind11/operators.h"

namespace py = pybind11;
using namespace pybind11::literals;

#define PYB11COMMA ,
#define PYB11LANGLE <
#define PYB11RANGLE >

""" % {"name" : name})

    # Includes
    if hasattr(modobj, "includes"):
        for inc in modobj.includes:
            ss('#include %s\n' % inc)
        ss("\n")

    # Use  namespaces
    if hasattr(modobj, "namespaces"):
        for ns in modobj.namespaces:
            ss("using namespace " + ns + ";\n")
        ss("\n")

    # Use objects from scopes
    if hasattr(modobj, "scopenames"):
        for scopename in modobj.scopenames:
            ss("using " + scopename + "\n")
        ss("\n")

    # Preamble
    if hasattr(modobj, "preamble"):
        if modobj.preamble:
            ss(modobj.preamble + "\n")
        ss("\n")

    # Some pybind11 types need their own preamble.
    for objname, obj in PYB11STLobjs(modobj):
        obj.preamble(modobj, ss, objname)
    ss("\n")

    # Trampolines
    PYB11generateModuleTrampolines(modobj, ss)

    # Declare the module
    ss("""
//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(%(name)s, m) {

""" % {"name"     : name,
      })

    if inspect.getdoc(modobj):
        ss("  m.doc() = ")
        for i, line in enumerate(inspect.getdoc(modobj).split('\n')):
            if i > 0:
                ss("            ")
            ss('"%s"\n' % line);
        ss("  ;\n")
    ss("\n")

    return

