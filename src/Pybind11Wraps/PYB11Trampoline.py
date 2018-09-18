#-------------------------------------------------------------------------------
# PYB11TrampolineGenerator
#-------------------------------------------------------------------------------
import inspect
import sys
import StringIO

from PYB11utils import *

#-------------------------------------------------------------------------------
# PYB11generateTrampoline
#
# Generate the trampoline class, including pure virtual hooks.
#-------------------------------------------------------------------------------
def PYB11generateTrampoline(klass, klassattrs, ssout):

    # Common preliminary code
    PYB11generateTrampolineClassStart(klass, klassattrs, ssout)

    # Bind the (unique) virtual methods for all classes up the inheritance tree.
    boundMethods = []
    for bklass in inspect.getmro(klass):
        bklassinst = bklass()
        bklassattrs = PYB11attrs(bklass)
        methods = [(mname, meth) for (mname, meth) in PYB11ClassMethods(bklass)
                   if not PYB11attrs(meth)["ignore"] and
                   (PYB11attrs(meth)["virtual"] or PYB11attrs(meth)["pure_virtual"])]
        for mname, meth in methods:
            
            # We build this method string up independent of the output stream
            # until we determine if it's already been generated.
            fs = StringIO.StringIO()
            ss = fs.write

            methattrs = PYB11attrs(meth)
            methattrs["returnType"] = eval("bklassinst." + mname + "()")
            assert methattrs["returnType"]    # We require the full spec for virtual methods
            ss("  virtual %(returnType)s %(cppname)s(" % methattrs)

            # Fill out the argument list for this method
            args = PYB11parseArgs(meth)
            for i, (argType, argName, default) in enumerate(args):
                ss(argType)
                if i < len(args) - 1:
                    ss(", ")
            if methattrs["const"]:
                ss(") const override { ")
            else:
                ss(") override { ")

            if methattrs["pure_virtual"]:
                ss("PYBIND11_OVERLOAD_PURE(%s, " % methattrs["returnType"])
            else:
                ss("PYBIND11_OVERLOAD(%s, " % methattrs["returnType"])
            ss(" %(namespace)s%(cppname)s," % klassattrs)
            if len(args) > 0:
                ss(" %(cppname)s, " % methattrs)
            else:
                ss(" %(cppname)s);" % methattrs)

            for i, (argType, argName, default) in enumerate(args):
                if i < nargs - 1:
                    ss(argName + ", ")
                else:
                    ss(argName + ");")
            ss(" }\n")

            # Is this a new binding?
            if not fs.getvalue() in boundMethods:
                ssout(fs.getvalue())
                boundMethods.append(fs.getvalue())
            fs.close()

    # Closing
    PYB11generateTrampolineClassEnd(klass, klassattrs, ssout)
    return

#-------------------------------------------------------------------------------
# PYB11generateTrampolineClassStart
#
# All the stuff up to the methods.
#-------------------------------------------------------------------------------
def PYB11generateTrampolineClassStart(klass, klassattrs, ssout):

    # Prepare in case there are templates lurking in here.
    fs = StringIO.StringIO()
    ss = fs.write

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Trampoline class for %(cppname)s
//------------------------------------------------------------------------------
#ifndef __trampoline_%(pyname)s__
#define __trampoline_%(pyname)s__

""" % klassattrs)

    # Namespaces
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("namespace " + ns + " {\n")

    # Class name
    ss("""
class PYB11Trampoline%(cppname)s: public %(cppname)s {
public:
  using %(cppname)s::%(cppname)s;   // inherit constructors

""" % klassattrs)

    # Finish up.
    ssout(fs.getvalue() % klassattrs["template_dict"])
    fs.close()

    return

#-------------------------------------------------------------------------------
# PYB11generateTrampolineClassEnd
#
# Finish up the code.
#-------------------------------------------------------------------------------
def PYB11generateTrampolineClassEnd(klass, klassattrs, ss):
    ss("};\n\n")
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("}\n")
    ss("\n#endif\n")
    return
