#-------------------------------------------------------------------------------
# PYB11Trampoline
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
def PYB11generateTrampoline(klass, ssout):

    # Check if this is a templated class definition
    klassattrs = PYB11attrs(klass)
    template_klass = len(klassattrs["template"]) > 0
    klassattrs["full_cppname"], klassattrs["mangle_cppname"] = klassattrs["cppname"], klassattrs["cppname"]
    if template_klass:
        stuff = PYB11mangle(klassattrs["template"])
        klassattrs["full_cppname"] += stuff[0]
        klassattrs["mangle_cppname"] += stuff[1]

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

    if template_klass:
        ss("template<")
        for i, name in enumerate(klassattrs["template"]):
            if i < len(klassattrs["template"]) - 1:
                ss("typename %s, " % name)
            else:
                ss("typename %s>\n" % name)

    # Class name
    ss("""class PYB11Trampoline%(cppname)s: public %(full_cppname)s {
public:
  using %(full_cppname)s::%(cppname)s;   // inherit constructors

""" % klassattrs)

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
            fms = StringIO.StringIO()
            ms = fms.write

            methattrs = PYB11attrs(meth)
            methattrs["returnType"] = eval("bklassinst." + mname + "()")
            assert methattrs["returnType"]    # We require the full spec for virtual methods
            ms("  virtual %(returnType)s %(cppname)s(" % methattrs)

            # Fill out the argument list for this method
            args = PYB11parseArgs(meth)
            for i, (argType, argName, default) in enumerate(args):
                ms(argType)
                if i < len(args) - 1:
                    ms(", ")
            if methattrs["const"]:
                ms(") const override { ")
            else:
                ms(") override { ")

            if methattrs["pure_virtual"]:
                ms("PYBIND11_OVERLOAD_PURE(%s, " % methattrs["returnType"])
            else:
                ms("PYBIND11_OVERLOAD(%s, " % methattrs["returnType"])
            ms(" %(namespace)s%(full_cppname)s," % klassattrs)
            if len(args) > 0:
                ms(" %(cppname)s, " % methattrs)
            else:
                ms(" %(cppname)s);" % methattrs)

            for i, (argType, argName, default) in enumerate(args):
                if i < nargs - 1:
                    ms(argName + ", ")
                else:
                    ms(argName + ");")
            ms(" }\n")

            # Is this a new binding?
            if not fms.getvalue() in boundMethods:
                ss(fms.getvalue())
                boundMethods.append(fms.getvalue())
            fms.close()

    # Closing
    ss("};\n\n")
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("}\n")
    ss("\n#endif\n")

    # Sub any template parameters.
    Tdict = {}
    for p in klassattrs["template"]:
        Tdict[p] = p
    ssout(fs.getvalue() % Tdict)

    return

