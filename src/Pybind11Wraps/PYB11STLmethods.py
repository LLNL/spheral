#-------------------------------------------------------------------------------
# PYB11STLmethods
#
# Thin wrappers to generate the pybind11 STL functions.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# std::vector
#-------------------------------------------------------------------------------
class PYB11_bind_vector:

    def __init__(self,
                 element,            # template element type for vector
                 opaque = False):    # should we make the type opaque
        self.element = element
        self.opaque = opaque
        return

    def preamble(self, modobj, ss, name):
        if self.opaque:
            ss('PYBIND11_MAKE_OPAQUE(std::vector<' + self.element + '>);\n')
        return

    def __call__(self, modobj, ss, name):
        ss('py::bind_vector<std::vector<' + self.element + '>>(m, "' + name + '");\n')
        return

#-------------------------------------------------------------------------------
# std::map
#-------------------------------------------------------------------------------
class PYB11_bind_map:

    def __init__(self,
                 key,                # template key type
                 value,              # template value type
                 opaque = False):    # should we make the container opaque
        self.key = key
        self.value = value
        self.opaque = opaque
        return

    def preamble(self, modobj, ss, name):
        if self.opaque:
            mangle_key = self.key.replace("::", "_")
            mangle_val = self.value.replace("::", "_")
            ttype = 'map_' + mangle_key + '_' + mangle_val
            ss('''
typedef std::map<%(key)s, %(value)s> %(ttype)s;
PYBIND11_MAKE_OPAQUE(%(ttype)s);
''' % {"key": self.key, "value": self.value, "ttype": ttype})
        return

    def __call__(self, modobj, ss, name):
        ss('py::bind_map<std::map<' + self.key + ', ' + self.value + '>>(m, "' + name + '");\n')
        return
