#-------------------------------------------------------------------------------
# PYB11Trampoline
#-------------------------------------------------------------------------------
import inspect
import sys
import StringIO

from PYB11utils import *

#-------------------------------------------------------------------------------
# PYB11generateModuleTrampolines
#
# Generate trampolines for any classes with virtual methods.
#-------------------------------------------------------------------------------
def PYB11generateModuleTrampolines(modobj, ss):
    klasses = PYB11classes(modobj)

    # Cull to just classes with virtual methods.
    klasses = [(name, klass) for (name, klass) in klasses if PYB11virtualClass(klass)]

    # Cull for things we're ignoring
    newklasses = []
    known_trampolines = []
    for name, klass in klasses:
        klassattrs = PYB11attrs(klass)
        template_klass = len(klassattrs["template"]) > 0
        mods = klassattrs["module"]
        if ((template_klass or not klassattrs["ignore"]) and                 # ignore flag (except for template class)?
            (klassattrs["pyname"] not in known_trampolines) and              # has this trampoline been generated?
            ((klass not in mods) or mods[klass] == modobj.PYB11modulename)): # is this class imported from another mod?
            newklasses.append((name, klass))
            known_trampolines.append(klassattrs["pyname"])
    klasses = newklasses

    # Generate trampolines
    for kname, klass in klasses:
        PYB11generateTrampoline(klass, ss)
    return

#-------------------------------------------------------------------------------
# PYB11generateTrampoline
#
# Generate the trampoline class, including pure virtual hooks.
#-------------------------------------------------------------------------------
def PYB11generateTrampoline(klass, ssout):

    klassattrs = PYB11attrs(klass)
    template_klass = len(klassattrs["template"]) > 0
    # if klassattrs["ignore"] and not template_klass:
    #     return

    # # This is a bit of trickery to let us use inheritance without regenerating trampolines
    # if"__known_trampolines" not in PYB11generateTrampoline.__dict__:
    #     PYB11generateTrampoline.__known_trampolines = []
    # if klassattrs["pyname"] in PYB11generateTrampoline.__known_trampolines:
    #     return
    # PYB11generateTrampoline.__known_trampolines.append(klassattrs["pyname"])

    # # We also screen out classes imported from other modules.
    # mods = klassattrs["module"]
    # if klass in mods:
    #     return

    # Prepare in case there are templates lurking in here.
    fs = StringIO.StringIO()
    ss = fs.write

    # Build the dictionary of template substitutions.
    Tdict = {key:key for key in klassattrs["template"]}
    if klassattrs["template_dict"]:
        for key in klassattrs["template_dict"]:
            if not key in Tdict:
                Tdict[key] = klassattrs["template_dict"][key]

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
                if name in klassattrs["template"]:
                    bklassname += name
                else:
                    if not name in klassattrs["template_dict"]:
                        raise RuntimeError, "Trampoline template base class error: %s is missing from specified template parameters %s\n  (class, base) = (%s, %s)" % (name, klassattrs["template_dict"], klass, bklass)
                    bklassname += klassattrs["template_dict"][name]
                if i < len(bklassattrs["template"]) - 1:
                    bklassname += ", "
            bklassname += ">"
        bklassnames.append(bklassname)
    assert len(bklassnames) == len(inspect.getmro(klass))
    bklassnames[0] = "PYB11self"

    # Class name
    ss("""class PYB11Trampoline%(cppname)s: public %(full_cppname)s {
public:
  using %(full_cppname)s::%(cppname)s;   // inherit constructors
  typedef %(full_cppname)s PYB11self;    // Necessary to protect macros below from names with commas in them
""" % klassattrs)
    for bklassname in bklassnames[1:]:
        if bklassname != PYB11mangle(bklassname):
            ss("  typedef %s %s;\n" % (bklassname, PYB11mangle(bklassname)))

    # # Use any nested class definitions
    # klasses = [(x, eval("klass.%s" % x)) for x in dir(klass) if (inspect.isclass(eval("klass.%s" % x)) and x in klass.__dict__)]
    # for (kname, nklass) in klasses:
    #     nklassattrs = PYB11attrs(nklass)
    #     ss("  typedef typename %(full_cppname)s::" % klassattrs)
    #     ss("%(cppname)s %(cppname)s;\n" % nklassattrs)

    # Any typedefs?
    if hasattr(klass, "PYB11typedefs"):
        typedefs = str(klass.PYB11typedefs)
    else:
        typedefs = ""

    # Bind the (unique) virtual methods for all classes up the inheritance tree.
    # We use an independent StringIO object for this, since we may have some new typedefs that
    # need to be added before this stuff is output to the source.
    methfms = StringIO.StringIO()
    boundMethods = []
    for (bklass, bklassname) in zip(inspect.getmro(klass), bklassnames):

        bklassinst = bklass()
        bklassattrs = PYB11attrs(bklass)
        methods = [(mname, meth) for (mname, meth) in PYB11ClassMethods(bklass)
                   if (not PYB11attrs(meth)["ignore"] and
                       (PYB11attrs(meth)["virtual"] or PYB11attrs(meth)["pure_virtual"]) and
                       mname in bklass.__dict__)]

        # Look for any template parameters of the base not shared by the class in question
        bklasssubs = {}
        for name in bklassattrs["template"]:
            if not name in klassattrs["template"]:
                assert name in klassattrs["template_dict"]
                bklasssubs[name] = klassattrs["template_dict"][name]

        for mname, meth in methods:
            
            # We build this method string up independent of the output stream
            # until we determine if it's already been generated.
            fms = StringIO.StringIO()

            methattrs = PYB11attrs(meth)
            methattrs["returnType"] = eval("bklassinst." + mname + "()")
            assert methattrs["returnType"]    # We require the full spec for virtual methods
            fms.write("  virtual %(returnType)s %(cppname)s(" % methattrs)

            # Fill out the argument list for this method
            args = PYB11parseArgs(meth)
            for i, (argType, argName, default) in enumerate(args):
                fms.write("%s %s" % (argType, argName))
                if i < len(args) - 1:
                    fms.write(", ")
            if methattrs["const"]:
                fms.write(") const override { ")
            else:
                fms.write(") override { ")

            # At this point we can make the call of whether this is a new method.
            try:
                thpt = fms.getvalue() % Tdict
            except:
                raise RuntimeError, "Unable to generate call descriptor for %s in %s->%s" % (mname, str(klass), bklassname)
            if not thpt in boundMethods:
                boundMethods.append(fms.getvalue() % Tdict)

                # Check if the returnType C++ name will choke PYBIND11_OVERLOAD*
                returnType = methattrs["returnType"]
                if PYB11badchars(returnType):
                    returnType = PYB11mangle(returnType)
                    typedefstring = "    typedef %s %s;\n" % (methattrs["returnType"], returnType)
                    if typedefstring not in typedefs:
                        typedefs += typedefstring
                    methattrs["returnType"] = returnType

                if methattrs["pure_virtual"]:
                    fms.write("PYBIND11_OVERLOAD_PURE(%(returnType)s, PYB11self, %(cppname)s, " % methattrs)
                else:
                    # HACK!  To workaround what appears to be a bug in overloading virtual method callbacks
                    # in pybind11 (see https://github.com/pybind/pybind11/issues/1547), we have to give
                    # the address of the object that actually implements it.  This is clealy not how a human
                    # should have to handle this, but since we're code generating this we can do this explicit
                    # workaround.
                    #fms.write("PYBIND11_OVERLOAD(%(returnType)s, PYB11self, %(cppname)s, " % methattrs)
                    fms.write("PYBIND11_OVERLOAD(%(returnType)s, " % methattrs)
                    fms.write(PYB11mangle(bklassname) + ", ")
                    fms.write("%(cppname)s, " % methattrs)

                for i, (argType, argName, default) in enumerate(args):
                    if i < len(args) - 1:
                        fms.write(argName + ", ")
                    else:
                        fms.write(argName)
                fms.write("); }\n")

                # Write to the method overloading stream.
                methfms.write(fms.getvalue())
            fms.close()

    # Write the full typdefs
    ss(typedefs + "\n")

    # Write the method overloads
    ss(methfms.getvalue())
    methfms.close()

    # Closing
    ss("};\n\n")
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("}\n")
    ss("\n#endif\n")

    # Sub any template parameters.
    ssout(fs.getvalue() % Tdict)

    return
