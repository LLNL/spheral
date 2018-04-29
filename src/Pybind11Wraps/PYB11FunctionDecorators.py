#-------------------------------------------------------------------------------
# Decorators for functions.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator as PYB11decorator           # To preserve wrapped functions args

#-------------------------------------------------------------------------------
# cppname
#-------------------------------------------------------------------------------
def PYB11cppname(f, arg):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11cppname = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# pyname
#-------------------------------------------------------------------------------
def PYB11cppname(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11pyname = True
    return PYB11decorator.decorate(f, wrapper)

