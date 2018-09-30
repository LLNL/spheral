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

    klassattrs = PYB11attrs(klass)
    template_klass = len(klassattrs["template"]) > 0
    if klassattrs["ignore"] and not template_klass:
        return

    # This is a bit of trickery to let us use inheritance without regenerating trampolines
    if"__known_trampolines" not in PYB11generateTrampoline.__dict__:
        PYB11generateTrampoline.__known_trampolines = []
    if klassattrs["pyname"] in PYB11generateTrampoline.__known_trampolines:
        return
    PYB11generateTrampoline.__known_trampolines.append(klassattrs["pyname"])

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

    # Build the base class hierarchy names
    bklassnames = []
    for bklass in inspect.getmro(klass):
        bklassattrs = PYB11attrs(bklass)
        bklassname = "%(namespace)s%(cppname)s" % bklassattrs
        if len(bklassattrs["template"]) > 0:
            bklassname += "<"
            for i, name in enumerate(bklassattrs["template"]):
                bklassname += name
                if i < len(bklassattrs["template"]) - 1:
                    bklassname += ", "
            bklassname += ">"
        bklassnames.append(bklassname)
    assert len(bklassnames) == len(inspect.getmro(klass))

    # Class name
    ss("""class PYB11Trampoline%(cppname)s: public %(full_cppname)s {
public:
  using %(full_cppname)s::%(cppname)s;   // inherit constructors
  typedef %(full_cppname)s PYB11self;    // Necessary to protect macros below from names with commas in them
""" % klassattrs)
    for bklassname in bklassnames:
        ss("  typedef %s %s;\n" % (bklassname, PYB11mangle(bklassname)))

    # Any typedefs?
    if hasattr(klass, "typedefs"):
        ss(klass.typedefs + "\n")

    # Bind the (unique) virtual methods for all classes up the inheritance tree.
    boundMethods = []
    for (bklass, bklassname) in zip(inspect.getmro(klass), bklassnames):

        bklassinst = bklass()
        bklassattrs = PYB11attrs(bklass)
        methods = [(mname, meth) for (mname, meth) in PYB11ClassMethods(bklass)
                   if (not PYB11attrs(meth)["ignore"] and
                       (PYB11attrs(meth)["virtual"] or PYB11attrs(meth)["pure_virtual"]) and
                       mname in bklass.__dict__)]

        # # HACK!
        # bklassname = "%(namespace)s%(cppname)s" % bklassattrs
        # if len(bklassattrs["template"]) > 0:
        #     bklassname += "<"
        #     for i, name in enumerate(bklassattrs["template"]):
        #         bklassname += name
        #         if i < len(bklassattrs["template"]) - 1:
        #             bklassname += ", "
        #     bklassname += ">"
        # bklassname = PYB11CPPsafe(bklassname)

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
                ms("%s %s" % (argType, argName))
                if i < len(args) - 1:
                    ms(", ")
            if methattrs["const"]:
                ms(") const override { ")
            else:
                ms(") override { ")

            # At this point we can make the call of whether this is a new method.
            if not fms.getvalue() in boundMethods:
                boundMethods.append(fms.getvalue())

                if methattrs["pure_virtual"]:
                    ms("PYBIND11_OVERLOAD_PURE(%(returnType)s, PYB11self, %(cppname)s, " % methattrs)
                else:
                    # HACK!  To workaround what appears to be a bug in overloading virtual method callbacks
                    # in pybind11 (see https://github.com/pybind/pybind11/issues/1547), we have to give
                    # the address of the object that actually implements it.  This is clealy not how a human
                    # should have to handle this, but since we're code generating this we can do this explicit
                    # workaround.
                    #ms("PYBIND11_OVERLOAD(%(returnType)s, PYB11self, %(cppname)s, " % methattrs)
                    ms("PYBIND11_OVERLOAD(%(returnType)s, " % methattrs)
                    ms(PYB11mangle(bklassname) + ", ")
                    ms("%(cppname)s, " % methattrs)

                for i, (argType, argName, default) in enumerate(args):
                    if i < len(args) - 1:
                        ms(argName + ", ")
                    else:
                        ms(argName)
                ms("); }\n")

                # Write to the out stream.
                ss(fms.getvalue())
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
