SPHERALSRC = /usr/workspace/wsrzd/owen/projects/Spheral/spheral/src
SPHERALBLD = /usr/workspace/wsrzd/owen/projects/Spheral/debug-pybind11-toss_3_x86_64_ib/Linux_x86_64
#SPHERALSRC = /home/owen/projects/Spheral/spheral/src
#SPHERALBLD = /home/owen/projects/Spheral/debug-pybind11/Linux_x86_64
PYTHON = $(SPHERALBLD)/bin/python
CXX = g++
CXXFLAGS = -fpic -fexceptions -std=c++11 -march=native -I.  -I$(SPHERALBLD)/include/python2.7  -I$(SPHERALBLD)/include

all:	testInheritance testTemplateInheritance

testInheritance:	testInheritance.py
	$(PYTHON) -c "import sys; sys.path += ['$(SPHERALSRC)/Pybind11Wraps']; from PYB11Generator import *; import testInheritance; PYB11generateModule(testInheritance)"
	$(CXX) $(CXXFLAGS) -shared -o testInheritance.so testInheritance.cc

testTemplateInheritance:	testTemplateInheritance.py
	$(PYTHON) -c "import sys; sys.path += ['$(SPHERALSRC)/Pybind11Wraps']; from PYB11Generator import *; import testTemplateInheritance; PYB11generateModule(testTemplateInheritance)"
	$(CXX) $(CXXFLAGS) -shared -o testTemplateInheritance.so testTemplateInheritance.cc

clean:
	rm -fr *.cc *.so
