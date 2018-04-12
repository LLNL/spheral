#-------------------------------------------------------------------------------
# Decorators for attributes to classes.
#-------------------------------------------------------------------------------
from functools import wraps

#-------------------------------------------------------------------------------
# Singleton
#-------------------------------------------------------------------------------
class singleton:

    def __init__(self, state=True):
        self.state = state
        return

    def __call__(self, cls):
        @wraps(cls)
        class Wrapped:
            cls.__pyb11_singleton = self.state
        return Wrapped

#-------------------------------------------------------------------------------
# Virtual method
#-------------------------------------------------------------------------------
def virtual(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        f.__pyb11_virtual = True
    return wrapper

#-------------------------------------------------------------------------------
# Pure virtual method
#-------------------------------------------------------------------------------
def pure_virtual(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        f.__pyb11_pure_virtual = True
    return wrapper

#-------------------------------------------------------------------------------
# const method
#-------------------------------------------------------------------------------
def const(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        f.__pyb11_const = True
    return wrapper

