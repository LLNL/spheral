#-------------------------------------------------------------------------------
# PYB11Publicist
#-------------------------------------------------------------------------------
import inspect
import sys
import StringIO

from PYB11utils import *

#-------------------------------------------------------------------------------
# PYB11generateModulePublicists
#
# Generate publicists for any classes with protected methods.
#-------------------------------------------------------------------------------
def PYB11generateModulePublicists(modobj, ss):
    klasses = PYB11classes(modobj)
    known_publicists = []
    for kname, klass in klasses:
        klassattrs = PYB11attrs(klass)
        template_klass = len(klassattrs["template"]) > 0
        mods = klassattrs["module"]
        if (PYB11protectedClass(klass) and
            ((template_klass or not klassattrs["ignore"]) and                 # ignore flag (except for template class)?
             (klassattrs["pyname"] not in known_publicists) and               # has this trampoline been generated?
            ((klass not in mods) or mods[klass] == modobj.PYB11modulename))): # is this class imported from another mod?
            PYB11generatePublicist(klass, ss)
    return

#-------------------------------------------------------------------------------
# PYB11generatePublicist
#
# Generate the publicist class, including pure virtual hooks.
#-------------------------------------------------------------------------------
def PYB11generatePublicist(klass, ssout):

    klassattrs = PYB11attrs(klass)
    template_klass = len(klassattrs["template"]) > 0
    if klassattrs["ignore"] and not template_klass:
        return

    # Prepare in case there are templates lurking in here.
    fs = StringIO.StringIO()
    ss = fs.write

    # Build the dictionary of template substitutions.
    Tdict = {}
    for p in klassattrs["template"]:
        Tdict[p] = p

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Publicist class for %(cppname)s
//------------------------------------------------------------------------------
#ifndef __publicist_%(pyname)s__
#define __publicist_%(pyname)s__

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
    ss("""class PYB11Publicist%(cppname)s: public %(full_cppname)s {
public:
""" % klassattrs)

    # Any typedefs?
    if hasattr(klass, "PYB11typedefs"):
        ss(klass.PYB11typedefs + "\n")

    # Publish the virtual methods of this class.
    methods = [(mname, meth) for (mname, meth) in PYB11ClassMethods(klass)
               if (PYB11attrs(meth)["protected"] and mname in klass.__dict__)]

    boundmeths = []
    for mname, meth in methods:
        methattrs = PYB11attrs(meth)
        if methattrs["cppname"] not in boundmeths:
            boundmeths.append(methattrs["cppname"])
            ss("  using %(full_cppname)s::" % klassattrs)
            ss("%(cppname)s;\n" % methattrs)

    # Closing
    ss("};\n\n")
    for ns in klassattrs["namespace"].split("::")[:-1]:
        ss("}\n")
    ss("\n#endif\n")

    # Sub any template parameters.
    ssout(fs.getvalue() % Tdict)

    return
