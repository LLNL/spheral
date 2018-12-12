from PYB11Generator import *

includes = ['<iostream>']

preamble = """
class Base {
public:
  Base()                             { printf("Base()\\n"); }
  virtual ~Base()                    {}
  virtual void funcBase() const      { printf("Base.funcBase()\\n"); }
  virtual void funcBaseOnly() const  { printf("Base.funcBaseOnly()\\n"); }
};

class A: public Base {
public:
  A()                                 { printf("A()\\n"); }
  virtual ~A()                        {}
  virtual void funcA(int x) const     { printf("A.funcA(%d)\\n", x); }
  virtual void funcAonly(int x) const { printf("A.funcAonly(%d)\\n", x); }
};

class B {
public:
  B()                                    { printf("B()\\n"); }
  virtual ~B()                           {}
  virtual void funcB(double x) const     { printf("B.funcB(%e)\\n", x); }
  virtual void funcBonly(double x) const { printf("B.funcBonly(%e)\\n", x); }
};

class C: public A, public B {
public:
  C()                                    { printf("C()"); }
  virtual ~C()                           {}
  virtual void funcA(int x) const        { printf("C.funcA(%d)\\n", -x); }
  virtual void funcB(double x) const     { printf("C.funcB(%e)\\n", -x); }
  virtual void funcConly(double x) const { printf("C.funcConly(%e)\\n", -x); }
};
"""

#-------------------------------------------------------------------------------
class Base:
    def pyinit(self):
        "Base default constructor"

    @PYB11virtual
    @PYB11const
    def funcBase(self):
        return "void"

    @PYB11virtual
    @PYB11const
    def funcBaseOnly(self):
        return "void"

#-------------------------------------------------------------------------------
class A(Base):
    def pyinit(self):
        "A default constructor"

    @PYB11virtual
    @PYB11const
    def funcA(self, x="int"):
        return "void"

    @PYB11virtual
    @PYB11const
    def funcAonly(self, x="int"):
        return "void"

#-------------------------------------------------------------------------------
class B:
    def pyinit(self):
        "A default constructor"

    @PYB11virtual
    @PYB11const
    def funcB(self, x="double"):
        return "void"

    @PYB11virtual
    @PYB11const
    def funcBonly(self, x="double"):
        return "void"

#-------------------------------------------------------------------------------
class C(A, B):
    def pyinit(self):
        "C default constructor"

    @PYB11virtual
    @PYB11const
    def funcA(self, x="int"):
        return "void"

    @PYB11virtual
    @PYB11const
    def funcB(self, x="double"):
        return "void"
