#-------------------------------------------------------------------------------
# Decorators for functions.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator as PYB11decorator           # To preserve wrapped functions args

#-------------------------------------------------------------------------------
# cppname
#-------------------------------------------------------------------------------
# def PYB11cppname(f, x):
#     def wrapper(f, *args, **kwargs):
#         return f(*args, **kwargs)
#     f.PYB11cppname = x
#     return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# pyname
#-------------------------------------------------------------------------------
def PYB11pyname(f, x):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11pyname = x
    return PYB11decorator.decorate(f, wrapper)

